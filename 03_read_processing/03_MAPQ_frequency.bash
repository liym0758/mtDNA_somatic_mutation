#!/bin/bash
#PBS -N MAPQ_Frequency
#PBS -j oe
#PBS -q batch
#PBS -o MAPQ_Frequency.out
#PBS -e MAPQ_Frequency.err
#PBS -l nodes=1:ppn=4

SPECIES=(
    #"human"
    "mouse"
    #"rat"
    #"chicken"
    #"cobra"
    #"sturgeon"
    #"stingray"
)

process_species() {
    local species=$1
    
    local work_dir="/PATH/TO/pbmc.mt/SNP/$species/mgatk.CleanBam/bam"
    cd "$work_dir" || exit
    
    local log_file="$work_dir/03_MAPQ_frequency.log"
    echo "===== $(date) ===== Start processing $species =====" > "$log_file"
    echo "Working directory: $work_dir" >> "$log_file"
    
    if [ ! -f "possorted_bam.chrM.bam" ]; then
        echo "ERROR: Input BAM file not found: possorted_bam.chrM.bam" | tee -a "$log_file"
        return 1
    fi
    
    if [ ! -f "possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.bam" ]; then
        echo "ERROR: Realigned BAM file not found: possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.bam" | tee -a "$log_file"
        return 1
    fi
    
    echo "Step 1: Extracting MAPQ values from original BAM..." | tee -a "$log_file"
    samtools view -@ 4 possorted_bam.chrM.bam 2>> "$log_file" | \
        awk '{print $5}' | sort | uniq -c | \
        awk '{print $2 "\t" $1}' > Original_mapq_frequency.txt 2>> "$log_file"
    
    if [ ! -s "Original_mapq_frequency.txt" ]; then
        echo "ERROR: Failed to create Original_mapq_frequency.txt" | tee -a "$log_file"
        return 1
    fi
    
    echo "Step 2: Extracting MAPQ values from realigned BAM..." | tee -a "$log_file"
    samtools view -@ 4 possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.bam 2>> "$log_file" | \
        awk '{print $5}' | sort | uniq -c | \
        awk '{print $2 "\t" $1}' > Final_mapq_frequency.txt 2>> "$log_file"
    
    if [ ! -s "Final_mapq_frequency.txt" ]; then
        echo "ERROR: Failed to create Final_mapq_frequency.txt" | tee -a "$log_file"
        return 1
    fi
    
    echo "Step 3: Generating MAPQ distribution plots with R..." | tee -a "$log_file"
    Rscript --vanilla -e '
    if(file.exists("Original_mapq_frequency.txt")) {
        data <- read.delim("Original_mapq_frequency.txt", header=FALSE, col.names=c("MAPQ", "Count"))
        png("Original_MAPQ_Distribution.png", width=1200, height=800, res=150)
        par(mar=c(5,5,4,2))
        barplot(data$Count, names.arg=data$MAPQ, 
                xlab="MAPQ Score", ylab="Number of Reads", 
                main="Original MAPQ Distribution",
                col=ifelse(data$MAPQ == 0, "red", 
                          ifelse(data$MAPQ == 60, "blue", "gray50")),
                border=NA, cex.names=0.8)
        abline(h=0, col="gray30")
        legend("topright", 
               legend=c("Unmapped (0)", "High Quality (60)", "Other"),
               fill=c("red", "blue", "gray50"), bty="n")
        box()
        dev.off()
        cat("Generated Original_MAPQ_Distribution.png\n")
    } else {
        cat("ERROR: Original_mapq_frequency.txt not found\n")
    }
    
    if(file.exists("Final_mapq_frequency.txt")) {
        data <- read.delim("Final_mapq_frequency.txt", header=FALSE, col.names=c("MAPQ", "Count"))
        png("Final_MAPQ_Distribution.png", width=1200, height=800, res=150)
        par(mar=c(5,5,4,2))
        barplot(data$Count, names.arg=data$MAPQ, 
                xlab="MAPQ Score", ylab="Number of Reads", 
                main="Final MAPQ Distribution",
                col=ifelse(data$MAPQ == 0, "red", 
                          ifelse(data$MAPQ == 60, "blue", "gray50")),
                border=NA, cex.names=0.8)
        abline(h=0, col="gray30")
        legend("topright", 
               legend=c("Unmapped (0)", "High Quality (60)", "Other"),
               fill=c("red", "blue", "gray50"), bty="n")
        box()
        dev.off()
        cat("Generated Final_MAPQ_Distribution.png\n")
    } else {
        cat("ERROR: Final_mapq_frequency.txt not found\n")
    }
    ' >> "$log_file" 2>&1
    
    if [ ! -f "Original_MAPQ_Distribution.png" ] || [ ! -f "Final_MAPQ_Distribution.png" ]; then
        echo "WARNING: One or more plot files not created" | tee -a "$log_file"
    else
        echo "Successfully generated MAPQ distribution plots" | tee -a "$log_file"
    fi
    
    echo "===== $(date) ===== Finished processing $species =====" | tee -a "$log_file"
    
    return 0
}

for species in "${SPECIES[@]}"; do
    start_time=$(date +%s)
    echo "Starting processing for $species at $(date)"
    
    process_species "$species"
    result=$?
    
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    
    if [ $result -eq 0 ]; then
        echo "Successfully processed $species in $runtime seconds"
    else
        echo "Failed to process $species after $runtime seconds. Check log: /PATH/TO/pbmc.mt/SNP/$species/mgatk.new/bam/03_MAPQ_frequency.log"
    fi
done

echo "All species processed!"
