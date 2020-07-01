import os
import subprocess as sp
import glob
from collections import Counter

camisim_path       = "/Users/gabrielraulet/Development/SummerWork/CAMISIM/metagenomesimulation.py"
pbsim_path         = "/Users/gabrielraulet/Development/SummerWork/CAMISIM/tools/PBSIM-PacBio-Simulator/src/pbsim"
error_profile_path = "/Users/gabrielraulet/Development/SummerWork/CAMISIM/tools/PBSIM-PacBio-Simulator/data"
samtools_path      = "/Users/gabrielraulet/Development/SummerWork/CAMISIM/tools/samtools-1.3/samtools"
temp_path          = "/Users/gabrielraulet/Development/SummerWork/metagenome_pipeline/tmp"
ncbi_taxdump_path  = "/Users/gabrielraulet/Development/SummerWork/metagenome_pipeline/input_configurations/taxonomy.tar.gz"
strainsim_path     = "/Users/gabrielraulet/Development/SummerWork/CAMISIM/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir"

def generate_config(output_directory, metadata_file, id_to_genome_file, size, fragments_size_mean, fragments_size_sd, num_samples):
    main_str             = "[Main]\nseed=632741178\nphase=0\nmax_processors=4\ndataset_id=DATASET\noutput_directory={}\ntemp_directory={}\ngsa=True\npooled_gsa=True\nanonymous=False\ncompress=1\n".format(output_directory, temp_path)
    readsim_str          = "\n[ReadSimulator]\nreadsim={}\nerror_profiles={}\nsamtools={}\nprofile=\nsize={}\ntype=pbsim\nfragments_size_mean={}\nfragments_size_standard_deviation={}\n".format(pbsim_path, error_profile_path, samtools_path, size, fragments_size_mean, fragments_size_sd)
    community_design_str = "\n[CommunityDesign]\ndistribution_file_paths=\nncbi_taxdump={}\nstrain_simulation_template={}\nnumber_of_samples={}\n".format(ncbi_taxdump_path, strainsim_path, num_samples)
    community0_str       = "\n[community0]\nmetadata={}\nid_to_genome_file={}\nid_to_gff_file=\ngenomes_total=25\ngenomes_real=25\nmax_strains_per_otu=1\nratio=1\nmode=replicates\nlog_mu=1\nlog_sigma=1\ngauss_mu=1\ngauss_sigma=1\nview=False\n".format(metadata_file, id_to_genome_file)

    with open("test_config.ini", "w") as config:
        config.write(main_str + readsim_str + community_design_str + community0_str)

def run_simulator():
    sp.call(["python", camisim_path, "test_config.ini"], shell=False)

def cat_reads(sample_folder, alignments_directory):
    sample_dir = sample_folder.rsplit("/", 1)[0]
    sample_folder_newname = sample_dir + "/sample_" + sample_folder.split("/")[-1].split("_sample_")[-1]
    os.rename(sample_folder, sample_folder_newname)
    sample_folder = sample_folder_newname
    with open(alignments_directory + "/reads/" + sample_folder.split("/")[-1] + "_reads.fq", "w") as reads_file:
        for zipped_fq in glob.glob(sample_folder + "/reads/*.fq.gz"):
            sp.call(["gunzip", zipped_fq], shell=False)
            with open(zipped_fq[:-3], "r") as zipped_file:
                for line in zipped_file.readlines():
                    if line.startswith("@"):
                        reads_file.write(line.rstrip() + "_" + sample_folder.split("/")[-1].split("_sample_")[-1] + '\n')
                    else:
                        reads_file.write(line)

def generate_overlap_alignments(alignments_directory):
    reads_directory = alignments_directory + "/reads"
    overlaps_directory = alignments_directory + "/overlaps"
    reads_files = [reads_file for reads_file in os.listdir(reads_directory) if reads_file.endswith(".fq")]
    for reads_file1 in reads_files:
        for reads_file2 in reads_files:
            reads_file1_abspath = reads_directory + "/" + reads_file1
            reads_file2_abspath = reads_directory + "/" + reads_file2
            overlap_path = overlaps_directory + "/" + reads_file1.split(".fq")[0] + "_vs_" + reads_file2.split(".fq")[0] + "_overlap.paf"
            sp.call(["minimap2", "-x", "ava-pb", reads_file1_abspath, reads_file2_abspath, "-o", overlap_path])

def paf_to_network(path_to_paf, network):
    with open(path_to_paf, "r") as paf:
        for line in paf.readlines():
            line_items = line.rstrip().split("\t")
            query_name = line_items[0]
            target_name = line_items[5]
            divergence_score = float(line.split("dv:f:")[1].split("\t")[0])
            network.write("{}\t{}\t{}\n".format(query_name, target_name, 1 - divergence_score))

def hipmcl_cluster(working_directory):
    sp.call(["docker", "run", "-it", "-d", "--name", "hipmcl_run", "--mount", "type=bind,source={}/docker_volume_folder,target=/data".format(working_directory),
             "hipmcl_image"])
    sp.call(["docker", "exec", "-it", "hipmcl_run", "sh", "-c", "\"/data/cluster.sh\""])
    sp.call(["docker", "kill", "hipmcl_run"])

def parse_cluster(cluster_list):
    counter_output = open("run.txt", "a")
    genomes = []
    for entry in cluster_list:
        genome_name = entry.split("-")[0]
        sample_name = entry.split("_", 2)[-1]
        genomes.append(genome_name)
    counter_output.write(str(Counter(genomes)) + "\n")
    counter_output.close()


def parse_clusters(network_filepath):
    with open(network_filepath, "r") as net:
        for line in net.readlines():
            words = line.rstrip().split()
            parse_cluster(words)



if __name__ == "__main__":
    working_directory = os.getcwd()
    #  output_directory  = working_directory + "/output"
    #  generate_config(output_directory,
                    #  os.getcwd() + "/input_configurations/arctic-species-metadata.tsv",
                    #  os.getcwd() + "/input_configurations/arctic-species.tsv",
                    #  0.7, 10000, 0, 2)

    #  run_simulator()

    #  alignments_directory = working_directory + "/alignments_output"
    #  os.mkdir(alignments_directory)
    #  os.mkdir(alignments_directory + "/reads")
    #  os.mkdir(alignments_directory + "/overlaps")
    #  sample_folders = glob.glob(output_directory + "/*_sample_*")
    #  for sample_folder in sample_folders:
        #  cat_reads(sample_folder, alignments_directory)
    #  generate_overlap_alignments(alignments_directory)

    #  with open("docker_volume_folder/network.out", "w") as network:
        #  for overlap_paf in glob.glob(alignments_directory + "/overlaps/*"):
            #  paf_to_network(overlap_paf, network)

    hipmcl_cluster(working_directory)
    parse_clusters(working_directory + "/docker_volume_folder/network.out.hipmcl")




