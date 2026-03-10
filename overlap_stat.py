import time
start_time = time.time()

from pybedtools import BedTool
import scipy.stats as stats

# Input files (change input file name)
g_quads = "pqsfinder_hg38_strong.bed"
histone_peaks = "diff5y_k27ac_ovr_H3_rep1_peaks.broadPeak"

# Output files (change output file name )
g4_histone_observed_overlaps = "observed_overlaps_pqsfinder_hg38_strong_WITH_diff5y_k27ac_ovr_H3_rep1_peaks.bed"

g4_intervals = BedTool(g_quads)
histone_intervals = BedTool(histone_peaks)

# Set the number of iterations
iterations = 1000

# Perform the randomized overlap analysis for multiple iterations
overlaps = []
for _ in range(iterations):
    random_intervals = histone_intervals.shuffle(genome='hg38', noOverlapping=True)
    overlap = len(g4_intervals.intersect(random_intervals, u=True))
    overlaps.append(overlap)

# Compare the observed number of overlaps with the randomized overlaps
observed_overlaps = g4_intervals.intersect(histone_intervals, u=True)
observed_overlap_count = len(observed_overlaps)
randomized_overlaps_above_observed = sum(overlap >= observed_overlap_count for overlap in overlaps)
average_overlap = sum(overlaps) / len(overlaps)
highest_overlap = max(overlaps)

# Save the observed overlaps to a BED file (change output file name )
observed_overlaps.saveas(g4_histone_observed_overlaps)

# Calculate standard deviation
std_dev = average_overlap / (iterations ** 0.5)

# Calculate z-score
z_score = (observed_overlap_count - average_overlap) / std_dev

# Inputs for hypergeometric test
genome_size = 3137300923 # https://www.ncbi.nlm.nih.gov/grc/human/data. Total non-N bases
peak_sizes = []
#total_size = 0
No_of_histone_peaks = len(histone_intervals)

# Calculate avg peak size
peak_sizes = [int(feature.fields[2]) - int(feature.fields[1]) for feature in g4_intervals]
average_peak_size = sum(peak_sizes) / len(peak_sizes)

## Arguments for hypergeometric test:
#total no of bins/peaks in the genome
total_pop = genome_size/average_peak_size
#no of histone_peaks in the genome
successful_in_total_pop = No_of_histone_peaks
#no of g_quads
sample_size = len(peak_sizes)
#no of observed_overlaps
successful_in_sample = observed_overlap_count

# hypergeometric test
p_value = stats.hypergeom.sf(successful_in_sample - 1, total_pop, successful_in_total_pop, sample_size)

# Print the results
print("no of g_quads:", sample_size)
print("no of histone_peaks in the genome:", successful_in_total_pop)
print("Observed Overlaps:", observed_overlap_count)
print("Iterations:", iterations)
print("Average Randomized Overlap:", average_overlap)
print("Randomized Overlaps above Observed:", randomized_overlaps_above_observed)
print("Highest Randomized Overlap:", highest_overlap)
print("Average size of each histone peak:", average_peak_size)
print("Total no of bins/peaks in the genome", total_pop)
print("z-score:", z_score)
print("p-value:", p_value)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time} seconds")
