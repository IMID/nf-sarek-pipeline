import os
import csv
import sys
from jinja2 import Environment, FileSystemLoader

def create_samplesheet(base_dir, output_dir, sample_params):
    seqtype, source, vendor, design, version, run = sample_params

    input_directory = os.path.join(base_dir, seqtype, source, vendor, design, version, run)
    output_directory = os.path.join(output_dir, seqtype, source, vendor, design, version, run)
    os.makedirs(output_directory, exist_ok=True)
    
    output_csv = os.path.join(output_directory, "samplesheet.csv")
    
    samplesheet_data = []
    
    for root, _, files in os.walk(input_directory):
        for file in files:
            if file.endswith(".fastq.gz"):
                filename = os.path.join(root, file)
                parts = file.split('_')
                
                if len(parts) >= 4:
                    patient = parts[0]
                    sample = parts[0]
                    lane = int(parts[2][1:])
                    read_type = parts[3]
                    
                    if read_type == "R1":
                        fastq_1 = filename
                        fastq_2 = filename.replace("R1", "R2")
                        
                        if os.path.exists(fastq_2):
                            samplesheet_data.append([patient, sample, lane, fastq_1, fastq_2])
    
    samplesheet_data.sort(key=lambda x: (x[0], x[2]))  # Sort by patient and lane
    
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["patient", "sample", "lane", "fastq_1", "fastq_2"])
        writer.writerows(samplesheet_data)
    
    print(f"Samplesheet CSV created: {output_csv}")

    return output_csv, output_directory

def create_yaml(samplesheet_file, output_dir, sample_params):
    seqtype, source, vendor, design, version, run = sample_params

    bed_dir = f"/ifs/data/SLURM/nf_pipelines/nf-sarek-pipeline/beds/{seqtype}/{source}/{vendor}/{design}/{version}"
    bed_file = os.path.join(bed_dir, f"{source}_{vendor}_{design}_{version}.bed")
    
    template_dir = "/ifs/data/SLURM/nf_pipelines/nf-sarek-pipeline/templates"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("alignment_variant_caling_template.yaml")
    
    yaml_content = template.render(
        samplesheet_file=samplesheet_file,
        output_dir=output_dir,
        bed_file=bed_file
    )
    
    # Define output YAML path
    yaml_output_path = os.path.join(output_dir, "params.yaml")
    
    # Write YAML file
    with open(yaml_output_path, "w") as yaml_file:
        yaml_file.write(yaml_content)
    
    print(f"YAML configuration file created: {yaml_output_path}")

def main():
    if len(sys.argv) != 7:
        print("Usage: python create_samplesheet.py <seq_type> <source> <vendor> <design> <version> <run>")
        sys.exit(1)
    
    base_dir = "/ifs/data/NGS/raw_data"
    output_dir = "/ifs/data/NGS/nf_sarek"
    sample_params = sys.argv[1:7]
    
    samplesheet_file, output_directory = create_samplesheet(base_dir, output_dir, sample_params)
    create_yaml(samplesheet_file, output_directory, sample_params)

if __name__ == "__main__":
    main()
