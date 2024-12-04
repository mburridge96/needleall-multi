use chrono::prelude::*;
use clap::Parser;
use indicatif::{ProgressBar, ProgressState, ProgressStyle};
use rocraters::ro_crate::constraints::*;
use rocraters::ro_crate::graph_vector::GraphVector;
use rocraters::ro_crate::{
    contextual_entity, data_entity, graph_vector, metadata_descriptor, rocrate, root,
};
use seq_io::fasta::{Reader, Record};
use std::collections::HashMap;
use std::env::current_exe;
use std::fmt::Write as FmtWrite;
use std::fs::OpenOptions;
use std::fs::{create_dir_all, File};
use std::io::{self, BufRead, BufReader, Write as IoWrite};
use std::path::PathBuf;
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use tempfile::NamedTempFile;

/// Cli tool for multi-process Emboss NeedleAll
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to target fasta file
    #[arg(short, long)]
    fasta: String,

    /// Working directory. Default value creates a wd in %Y%m%d%H%M%S format
    #[arg(short, long, default_value = "date/")]
    working_dir: String,

    /// Identities output file name
    #[arg(short, long, default_value = "identities.tsv")]
    outfile: String,

    /// Needle all error file name
    #[arg(short, long, default_value = "needle_error.error")]
    errorfile: String,

    /// Gap open penalty
    #[arg(short, long, default_value = "10.0")]
    gap_open_penalty: f32,

    /// Gap extend penalty
    #[arg(short, long, default_value = "0.5")]
    gap_extend_penalty: f32,

    /// Threshold for result
    #[arg(short, long, default_value = "-1.0")]
    threshold: f32,

    /// Set cpu number for multi-processing
    #[arg(short, long, default_value = "0")]
    cpu_count: usize,

    /// Time debug
    #[arg(short, long, default_value = "false")]
    debug_time: bool,

    /// RO-Crate
    #[arg(short, long, default_value = "true")]
    rocrate: bool,

    // Agent
    #[arg(short, long, default_value = "Unknown")]
    agent: String,
}

fn main() {
    let args = Args::parse();

    let global_run_time = Instant::now();
    let start_time = Utc::now();

    // Check wd
    let wd: PathBuf = if args.working_dir == "date/" {
        let wd_name: String = Local::now().format("%Y%m%d%H%M%S").to_string();
        PathBuf::from(wd_name)
    } else {
        PathBuf::from(&args.working_dir)
    };

    // Make sure the working directory is made
    if let Err(e) = create_dir_all(&wd) {
        eprintln!("Could not create directory {:?}: {}", wd, e);
    }

    // set cpu count
    if args.cpu_count != 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.cpu_count)
            .build_global()
            .unwrap();
    } else {
        let num_cpus = num_cpus::get();
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_cpus - 2)
            .build_global()
            .unwrap();
    }

    let needleall = NeedleAll {
        outfile: &wd.join(&args.outfile).display().to_string(),
        errorfile: &wd.join(&args.errorfile).display().to_string(),
        gap_open_penalty: args.gap_open_penalty,
        gap_extend_penalty: args.gap_extend_penalty,
        threshold: args.threshold,
        debug_time: args.debug_time,
        rocrate: args.rocrate,
        agent: &args.agent,
    };
    println!("Needleall input: {:?}", needleall);

    let seq_data = get_seq_data(&args.fasta);

    process_seq_data_parallel(&seq_data, &needleall);

    let end_time = Utc::now();
    let global_run_end = global_run_time.elapsed();
    println!("RUN TIME: {:.2}s", global_run_end.as_secs());
    build_crate(wd, &args, &needleall, start_time, end_time);
}

// process in parallel
fn process_seq_data_parallel(seq_data: &SeqData, needle: &NeedleAll) {
    use rayon::prelude::*;

    let total_reps = seq_data.reps.len();
    let file_lock = Arc::new(Mutex::new(()));
    let pb = ProgressBar::new(total_reps as u64);
    pb.set_style(ProgressStyle::with_template(
        "{spinner:.green} [{elapsed_precise}] {bar:50.cyan/blue} {pos:>7}/{len:7} {msg} ({eta})",
    )
    .unwrap()
    .with_key("eta", |state: &ProgressState, w: &mut dyn FmtWrite| write!(w, "{:.1}s", state.eta().as_secs_f64()).unwrap())
    .progress_chars("##-"));

    seq_data
        .reps
        .par_iter()
        .enumerate()
        .for_each(|(index, &_rep)| {
            let run_time = Instant::now();

            let job = make_job(index, total_reps);
            let files = make_fasta(seq_data, index, job[0] as usize);

            let needle_time = Instant::now();
            run_needle_all(&files, needle);
            let needle_time_end = needle_time.elapsed();

            let result = parse_needle_output(&files);

            let file_lock = Arc::clone(&file_lock);
            match result {
                Ok(values) => {
                    let _lock = file_lock.lock().unwrap();
                    write_identities(&values, &needle.outfile, &needle.threshold)
                }
                Err(e) => println!("Failed to write: {}", e),
            }
            let run_time_end = run_time.elapsed();
            let run_vs_needle = run_time_end - needle_time_end;
            if needle.debug_time {
                println!(
                    "R{} - Run time: {}, - Needall time: {} - File IO handicap {}",
                    _rep,
                    run_time_end.as_secs_f32(),
                    needle_time_end.as_secs_f32(),
                    run_vs_needle.as_secs_f32()
                )
            }
            pb.inc(1);
        });
}

#[derive(Debug)]
struct NeedleAll<'a> {
    outfile: &'a String,
    errorfile: &'a String,
    gap_open_penalty: f32,
    gap_extend_penalty: f32,
    threshold: f32,
    debug_time: bool,
    rocrate: bool,
    agent: &'a String,
}
#[derive(Debug)]
struct SeqData {
    reps: Vec<u32>,
    ids: Vec<Vec<u8>>,
    seqs: Vec<Vec<u8>>,
}

struct NeedleAllFiles {
    aseq: NamedTempFile,
    bseq: NamedTempFile,
    needle: NamedTempFile,
}

fn get_seq_data(path: &str) -> SeqData {
    let mut reader = Reader::from_path(path).unwrap();
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut reps: Vec<u32> = Vec::new();

    let mut count = 0;

    while let Some(record) = reader.next() {
        let record = record.expect("Error reading record");

        ids.push(record.head().to_vec());
        seqs.push(record.seq().to_vec());
        reps.push(count);
        count += 1;
    }

    SeqData { reps, ids, seqs }
}

fn make_job(rep: usize, total_reps: usize) -> Vec<u32> {
    let mut jobs: Vec<u32> = Vec::new();

    for i in rep..total_reps {
        jobs.push(i as u32);
    }
    jobs
}

fn run_needle_all(files: &NeedleAllFiles, needle: &NeedleAll) {
    let _output = Command::new("needleall")
        .arg("-asequence")
        .arg(files.aseq.path())
        .arg("-bsequence")
        .arg(files.bseq.path())
        .arg("-gapopen")
        .arg(needle.gap_open_penalty.to_string())
        .arg("-gapextend")
        .arg(needle.gap_extend_penalty.to_string())
        .arg("-outfile")
        .arg(files.needle.path())
        .arg("-aformat3")
        .arg("pair")
        .arg("-errfile")
        .arg(&needle.errorfile)
        .output()
        .expect("failed to execute process");
}

fn make_temp() -> NeedleAllFiles {
    let files = NeedleAllFiles {
        aseq: tempfile::Builder::new()
            .suffix(".fasta")
            .tempfile()
            .unwrap(),
        bseq: tempfile::Builder::new()
            .suffix(".fasta")
            .tempfile()
            .unwrap(),
        needle: tempfile::Builder::new()
            .suffix(".fasta")
            .tempfile()
            .unwrap(),
    };
    files
}

fn make_fasta(seq_data: &SeqData, index: usize, job: usize) -> NeedleAllFiles {
    let mut files = make_temp();
    let start_index: usize = job;

    files.aseq.write_all(b">").unwrap(); // ensure header line starts with '>'
    files.aseq.write_all(&seq_data.ids[index]).unwrap();
    files.aseq.write_all(b"\n").unwrap(); // newline after header
    files.aseq.write_all(&seq_data.seqs[index]).unwrap();
    files.aseq.write_all(b"\n").unwrap(); // newline after sequence

    for (id, seq) in seq_data
        .ids
        .iter()
        .skip(start_index)
        .zip(seq_data.seqs.iter().skip(start_index))
    {
        files.bseq.write_all(b">").unwrap(); // ensure header line starts with '>'
        files.bseq.write_all(id).unwrap();
        files.bseq.write_all(b"\n").unwrap(); // newline after header
        files.bseq.write_all(seq).unwrap();
        files.bseq.write_all(b"\n").unwrap(); // newline after sequence
    }
    files
}

/// Parse needle output
fn parse_needle_output(files: &NeedleAllFiles) -> io::Result<Vec<(String, String, f32)>> {
    let file = File::open(&files.needle)?;
    let reader = BufReader::new(file);

    let mut identities = Vec::new();
    let mut id1 = String::new();
    let mut id2 = String::new();
    let mut identity: f32;

    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("# 1:") {
            id1 = line.split_whitespace().nth(2).unwrap_or("").to_string();
        } else if line.starts_with("# 2:") {
            id2 = line.split_whitespace().nth(2).unwrap_or("").to_string();
        } else if line.starts_with("# Identity") {
            let parts: Vec<&str> = line
                .split_whitespace()
                .nth(2)
                .unwrap_or("")
                .split('/')
                .collect();
            if parts.len() == 2 {
                if let (Ok(num), Ok(den)) = (parts[0].parse::<f32>(), parts[1].parse::<f32>()) {
                    identity = num / den;
                    identities.push((id1.clone(), id2.clone(), identity));
                }
            }
        }
    }
    Ok(identities)
}

fn write_identities(identities: &[(String, String, f32)], outfile: &String, threshold: &f32) {
    let mut file = OpenOptions::new()
        .create(true)
        .append(true)
        .open(outfile)
        .expect("Unable to open file");

    for (id1, id2, identity) in identities {
        if identity >= threshold {
            writeln!(file, "{}\t{}\t{:.6}", id1, id2, identity).expect("Unable to write to file");
        }
    }
}

fn build_crate(
    wd: PathBuf,
    args: &Args,
    needle: &NeedleAll,
    start_time: DateTime<Utc>,
    end_time: DateTime<Utc>,
) {
    let mut rocrate = rocrate::RoCrate {
        context: rocrate::RoCrateContext::ReferenceContext(
            "https://w3id.org/ro/crate/1.1/context".to_string(),
        ),
        graph: Vec::new(),
    };

    // Set new MetadataDescriptor
    let description = metadata_descriptor::MetadataDescriptor {
        id: "ro-crate-metadata.json".to_string(),
        type_: DataType::Term("CreativeWork".to_string()),
        conforms_to: Id::Id("https://w3id.org/ro/crate/1.1".to_string()),
        about: Id::Id("./".to_string()),
        dynamic_entity: None,
    };

    rocrate
        .graph
        .push(GraphVector::MetadataDescriptor(description));

    // Create new RootDataEntity
    let root_data_entity = root::RootDataEntity {
        id: "./".to_string(),
        type_: DataType::Term("Dataset".to_string()),
        name: "Needleall".to_string(),
        description: "Needleall result for the generation of all-to-all sequence identities"
            .to_string(),
        date_published: Utc::now().to_string(),
        license: License::Id(Id::Id("MIT LICENSE".to_string())),
        dynamic_entity: None,
    };

    rocrate
        .graph
        .push(GraphVector::RootDataEntity(root_data_entity));

    rocrate
        .graph
        .push(GraphVector::DataEntity(data_entity::DataEntity {
            id: needle.outfile.to_string(),
            type_: DataType::Term("File".to_string()),
            dynamic_entity: Some(HashMap::from([
                (
                    "name".to_string(),
                    EntityValue::EntityString("Identities".to_string()),
                ),
                (
                    "description".to_string(),
                    EntityValue::EntityString(
                        "Needleman All Identities file for all-to-all similarity comparison"
                            .to_string(),
                    ),
                ),
            ])),
        }));

    rocrate
        .graph
        .push(GraphVector::DataEntity(data_entity::DataEntity {
            id: needle.errorfile.to_string(),
            type_: DataType::Term("File".to_string()),
            dynamic_entity: Some(HashMap::from([
                (
                    "name".to_string(),
                    EntityValue::EntityString("Needleman Error File".to_string()),
                ),
                (
                    "description".to_string(),
                    EntityValue::EntityString(
                        "Contains errors that occured during needleman".to_string(),
                    ),
                ),
            ])),
        }));

    rocrate.graph.push(GraphVector::ContextualEntity(
        contextual_entity::ContextualEntity {
            id: "#needleman_1".to_string(),
            type_: DataType::Term("CreateAction".to_string()),
            dynamic_entity: Some(HashMap::from([
                (
                    "name".to_string(),
                    EntityValue::EntityString("Needleall run".to_string()),
                ),
                (
                    "description".to_string(),
                    EntityValue::EntityString("Needleall run instance".to_string()),
                ),
                (
                    "agent".to_string(),
                    EntityValue::EntityId(Id::Id(needle.agent.to_string())),
                ),
                (
                    "startTime".to_string(),
                    EntityValue::EntityString(start_time.to_string()),
                ),
                (
                    "endTime".to_string(),
                    EntityValue::EntityString(end_time.to_string()),
                ),
                (
                    "instrument".to_string(),
                    EntityValue::EntityId(Id::Id(env!("CARGO_PKG_REPOSITORY").to_string())),
                ),
                (
                    "object".to_string(),
                    EntityValue::EntityId(Id::Id(args.fasta.to_string())),
                ),
                (
                    "result".to_string(),
                    EntityValue::EntityId(Id::IdArray(Vec::from([
                        needle.errorfile.to_string(),
                        needle.outfile.to_string(),
                    ]))),
                ),
            ])),
        },
    ));

    let entity_types = Vec::from(["File".to_string(), "SoftwareSourceCode".to_string()]);
    rocrate.graph.push(GraphVector::ContextualEntity(
        contextual_entity::ContextualEntity {
            id: env!("CARGO_PKG_REPOSITORY").to_string(),
            type_: DataType::TermArray(entity_types),
            dynamic_entity: Some(HashMap::from([
                (
                    "name".to_string(),
                    EntityValue::EntityString(env!("CARGO_PKG_NAME").to_string()),
                ),
                (
                    "description".to_string(),
                    EntityValue::EntityString(env!("CARGO_PKG_DESCRIPTION").to_string()),
                ),
                (
                    "version".to_string(),
                    EntityValue::EntityString(env!("CARGO_PKG_VERSION").to_string()),
                ),
                (
                    "codeRepository".to_string(),
                    EntityValue::EntityString(env!("CARGO_PKG_REPOSITORY").to_string()),
                ),
            ])),
        },
    ));

    rocraters::ro_crate::write::write_crate(
        &rocrate,
        wd.join("ro-crate-metadata.json").display().to_string(),
    )
}
