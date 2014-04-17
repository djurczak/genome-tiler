require 'command'
require_relative 'lib/genome_tiler.rb'

reference_file = "/clustertmp/brennecke/jurczak/projects/genome-tiles/data/dmel-all-chromosome-r5.51.noUex.fasta"
tokens = reference_file.split(".")
reference_file_rev = tokens[0..-2].join(".") + ".rev." + tokens[-1]

task :default do
  instance = GenomeTiler.new

  threads = []
  threads.push(Thread.new {
    fOut = File.open("/clustertmp/brennecke/jurczak/projects/genome-tiles/results/dmel-all-chromosome-r5.51.noUex.win32.fasta", 'w')
    instance.split_into_windows(reference_file, fOut, 32, {})
  })

  threads.push(Thread.new {
    fOut = File.open("/clustertmp/brennecke/jurczak/projects/genome-tiles/results/dmel-all-chromosome-r5.51.noUex.win20.fasta", 'w')
    instance.split_into_windows(reference_file, fOut, 20, {})
  })

  puts "Flipping reference"
  threads.push(Thread.new {
    File.open(reference_file_rev, 'w') do |fOut|
      instance.reverse_complement(reference_file, fOut)
    end

    t = Thread.new do |thread|
      fOut = File.open("/clustertmp/brennecke/jurczak/projects/genome-tiles/results/dmel-all-chromosome-r5.51.noUex.win20.rev.fasta", 'w')
      instance.split_into_windows(reference_file_rev, fOut, 32, {})
    end

    t2 = Thread.new do |thread|
      fOut = File.open("/clustertmp/brennecke/jurczak/projects/genome-tiles/results/dmel-all-chromosome-r5.51.noUex.win32.rev.fasta", 'w')
      instance.split_into_windows(reference_file_rev, fOut, 32, {})
    end

    t.join
    t2.join
  })

  threads.each(&:join)
  puts "Done."

  threads.clear

  puts "Generating BED files"
  threads.push(Thread.new {
    cmdline = ". /etc/profile.d/modules.sh && module load samtools bowtie && bowtie -f -v0 -a --best --strata --sam -p12 /groups/brennecke/annotation/genome/5_51/dmel-all-chromosome-r5.51.noUex /clustertmp/brennecke/jurczak/projects/genome-tiles/results/dmel-all-chromosome-r5.51.noUex.win32.rev.fasta | samtools view -Sb - | bamToBed -i - > dmel-all-chromosome-r5.51.noUex.win32.rev.bed"
    status = Command.run(cmdline)
    puts status.stderr unless status.success?
  })

  threads.push(Thread.new {
    cmdline = ". /etc/profile.d/modules.sh && module load samtools bowtie && bowtie -f -v0 -a --best --strata --sam -p12 /groups/brennecke/annotation/genome/5_51/dmel-all-chromosome-r5.51.noUex /clustertmp/brennecke/jurczak/projects/genome-tiles/results/dmel-all-chromosome-r5.51.noUex.win20.rev.fasta | samtools view -Sb - | bamToBed -i - > dmel-all-chromosome-r5.51.noUex.win20.rev.bed"
    status = Command.run(cmdline)
    puts status.stderr unless status.success?
  })

  threads.each(&:join)
  puts "done"
end
