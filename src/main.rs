use petgraph::stable_graph::StableGraph;
use petgraph::Directed;
use petgraph::graph::*;
use bio::io::fasta;
use std::collections::HashMap;
use petgraph::visit::EdgeRef;
use std::path::{Path, PathBuf};
use std::fs;
use anyhow::Error;
use std::ops::Deref;
use anyhow::anyhow;
use std::rc::Rc;
use bio::alphabets::dna::revcomp;

fn main() -> Result<(), Error> {
    let mut graph = StableGraph::<Vec<u8>, &str,Directed,usize>::default();

    let path = Path::new("/home/liam/ferraph/");
    let fns = get_fasta_path(&path)?;
    let mut seqs: Vec<Vec<Vec<usize>>> = Vec::new();
    let mut kmers: HashMap<Rc<Vec<u8>>, usize> = HashMap::new();
    for file in fns {
        let reader = fasta::Reader::from_file(&file)?;
        let mut gen_vec = Vec::new();
        for result in reader.records() {
           let mut seq_vec = Vec::new();
           let record = result.expect("Err");
           let overlap = record
                .seq()
                .windows(13)
                .map(|x| x.to_vec())
                .collect::<Vec<Vec<u8>>>();
            let mut isfirst = true;
            let mut lastkmer = NodeIndex::new(0);
            let mut lastrev = false;
            for km in overlap {
                let mut kmer = Rc::new(km);
                let rev = revcomp(Rc::deref(&kmer));
                let mut isrev = false;
                if Rc::deref(&kmer) > &rev {
                    kmer = Rc::new(rev);
                    isrev = true;
                }
                let edgetype;
                let nodeid;
                if !kmers.contains_key(&kmer) {
                    nodeid = graph.add_node(kmer.to_vec());
                    isfirst = false;
                    kmers.insert(Rc::clone(&kmer), nodeid.index());
                } else {
                    nodeid = NodeIndex::new(*kmers.get(&kmer).unwrap());
                }
                if !isfirst {
                    if lastrev == true {
                        if isrev == true {
                            edgetype = "RR";
                        } else {
                            edgetype = "RF";
                        }
                    } else {
                        if isrev == true {
                            edgetype = "FR";
                        } else {
                            edgetype = "FF";
                        }
                    }
                    let mut edgeid;
                    match &graph.find_edge(lastkmer, nodeid) {
                        Some(n) => {
                            edgeid = *n;
                            let mut newedge = true;
                            &graph.edges(lastkmer).for_each(|edge| {
                                if edge.target() == nodeid {
                                    if edge.weight() == &edgetype {
                                        edgeid = edge.id().clone();
                                        newedge = false;
                                    } else {
                                        println!("New Edge? {:#?} vs. {:#?}", edge.weight(), &edgetype)
                                    }
                                }
                            });

                            if newedge == false {
                                
                            } else {
                                println!("edgeid {:#?} from lastkmer {:#?} to nodeid {:#?}, oldedge {:#?}", &n, &lastkmer, nodeid, &graph[*n]);
                                edgeid = graph.add_edge(lastkmer, nodeid, edgetype);
                                println!("edgeid {:#?} from lastkmer {:#?} to nodeid {:#?}, newedge {:#?}", &edgeid, &lastkmer, &nodeid, &graph[edgeid]);
                            }
                        },
                        None => {
                            edgeid = graph.add_edge(lastkmer, nodeid, edgetype);
                        },
                    }
                    seq_vec.push(edgeid.index());
                }
                lastkmer = nodeid;
                lastrev = isrev;
            }
            gen_vec.push(seq_vec);
        }
        seqs.push(gen_vec);
    }
    println!("# Kmers: {:#?}", kmers.len());
    println!("# Nodes: {:#?}", graph.node_count());
    println!("# Edges {:#?}", graph.edge_count());
    println!("# Genomes: {:#?}", seqs[0].len());

    let checkind = EdgeIndex::new(811398);
    let connected = graph.edge_endpoints(checkind);
    match connected {
        Some(endnodes) => { print!("Seq {:#?} is connected to ", std::str::from_utf8(&graph[endnodes.0]));
                     print!("{:#?}", std::str::from_utf8(&graph[endnodes.1]))},
        None => println!("{:#?} ", "Error, edge is not connected to any nodes"),
    }
    println!("by a {:#?} edge", &graph[checkind]);

    Ok(())
}

fn get_fasta_path(file_or_dir: &Path) -> Result<Vec<PathBuf>, Error> {
    if file_or_dir.is_file() {
        Ok(vec![file_or_dir.to_path_buf()])
    } else if file_or_dir.is_dir() {
        let all_files = recurse_directory(file_or_dir)?;

        // Of all the files, we want to keep only the fasta files
        let mut fasta_files: Vec<PathBuf> = Vec::new();

        // To ensure we only get fasta files, we open each file, and attempt
        // to get the first fasta record of each. If we succeed, we add the file to
        // the vector of fasta files. If not, we do nothing.
        for f in all_files {
            let reader = fasta::Reader::from_file(&f);
            match reader {
                Ok(gf) => {
                    for record in gf.records() {
                        match record {
                            Ok(..) => fasta_files.push(f),
                            Err(..) => {}
                        }
                        break;
                    }
                }
                Err(..) => {}
            }
        }

        // Check to see if any fasta files were found. If not, return an error
        if fasta_files.is_empty() {
            Err(anyhow!(
                "No valid fasta files found"
            ))
        } else {
            Ok(fasta_files)
        }
    } else {
        Err(anyhow!("No valid files found"))
    }
}

// Path only holds a reference to the path string
// PathBuf owns the string
fn recurse_directory(p: &Path) -> Result<Vec<PathBuf>, Error> {
    let mut af: Vec<PathBuf> = Vec::new();
    for entry in fs::read_dir(p)? {
        let e = entry?;
        let path = e.path();

        if path.is_dir() {
            recurse_directory(&path)?;
        } else {
            af.push(path);
        }
    }
    Ok(af)
}