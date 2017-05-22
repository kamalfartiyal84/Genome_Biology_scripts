#!/usr/bin/env ruby

if ARGV.length < 3
  STDERR.puts "Usage: #{File.basename($0)} sample copy_number_file centromere_bin_file"
  exit
end

sample = ARGV.shift
copy_number_file = ARGV.shift
centromere_bin_file = ARGV.shift

centromeres = {}
centromere_start = {}
File.open(centromere_bin_file).each do |line|
	location = line.chomp.split(/\t/).first
	centromeres[location] = 1
	(chrom, start, stop) = location.split(/[:\-]/)
	centromere_start[chrom] = start.to_i unless centromere_start[chrom]
end

MAX_NA_COUNT = 50

current_chrom = nil
current_start = nil
current_arm = nil
current_stop = nil
current_segmented = nil
current_count = 0
na_count = 0

puts "sampleID	chrom	arm	start.pos	end.pos	n.probes	mean"

File.open(copy_number_file).each do |line|
	(location, chrom, start, stop, copynumber, segmented) = line.chomp.split(/\t/)
	next if location.eql?("location")
	start = start.to_i
	stop = stop.to_i

# do we need to write the previous segment?
#  - change of chromosome
#  - change of segmented copy number
#  - number of consecutive NA values above threshold
#  - in centromere

	if !chrom.eql?(current_chrom) or (!segmented.eql?("NA") and !segmented.eql?(current_segmented)) or na_count > MAX_NA_COUNT or centromeres.include?(location)

#		STDERR.puts "NA: #{current_chrom} #{current_stop + 1} #{start - 1} #{na_count}" if na_count > MAX_NA_COUNT

		puts [ sample, current_chrom, current_arm, current_start, current_stop, current_count, current_segmented ].join("\t") if current_chrom

		current_chrom = nil
		current_start = nil
		current_arm = nil
		current_stop = nil
		current_segmented = nil
		current_count = 0
		na_count = 0
	end

	if segmented.eql?("NA")
		na_count += 1
	else
		na_count = 0

		unless current_segmented
			current_chrom = chrom
			current_start = start
			current_arm = current_start < centromere_start[current_chrom] ? "p" : "q"
			current_segmented = segmented
			current_count = 0
		end

		current_stop = stop
		current_count += 1
	end

end

puts [ sample, current_chrom, current_arm, current_start, current_stop, current_count, current_segmented ].join("\t") if current_chrom

