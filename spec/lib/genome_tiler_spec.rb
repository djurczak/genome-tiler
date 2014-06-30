require_relative File.join('..', '..', 'lib', 'genome_tiler.rb')

describe GenomeTiler do
  let(:data) { File.join(File.dirname(__FILE__), '..', 'fixtures', 'sample_data.fasta') }

  before(:each) do
    @output = StringIO.new
  end

  describe "#reverse_complement" do
    let(:data) { File.join(File.dirname(__FILE__), '..', 'fixtures', 'two_sequences.fasta') }

    context "given a fasta file" do
      it "converts all sequences and stores them in the output file" do
        instance = GenomeTiler.new
        instance.reverse_complement(data, @output)
        expect(@output.string.split("\n").count).to eq 4
      end

      it "correctly flips all sequences" do
        instance = GenomeTiler.new

        instance.reverse_complement(data, @output)

        flipped = "ggaacatgaatgaaaataaaaacaagactggatcagcattacgtgaccct"
        @output.string.split("\n").each_with_index do |line,i|
          next if i.even?
          expect(line).to eq flipped
        end
      end
    end

    describe "#each_sequence" do
      it "yields control for every sequene in a fasta file" do
        instance = GenomeTiler.new
        expect { |b| instance.each_sequence(data, &b) }.to yield_control.exactly(2).times
      end
    end

    describe "#reverse_complement_sequence" do
      it "returns the reverse complement of a nucleotide sequence" do
        instance = GenomeTiler.new
        expect(instance.reverse_complement_sequence("AGCTTT").downcase).to eq "AAAGCT".downcase
      end
    end
  end


  describe "#split_into_windows" do
    it "takes a fasta file and splits its sequences into overlapping windows" do
      instance = GenomeTiler.new
      instance.split_into_windows(data, @output, 20, {})
      expect(@output.string.split("\n").count).to eq 62
    end

    describe "#each_window_in_data" do

      it "yields control for each possible split for a sequence" do
        instance = GenomeTiler.new
        expect {
          |b| instance.each_window_in_data(data, 20, &b)
        }.to yield_control.exactly(31).times
      end

      it "yields a position-specific subsequence" do
        i = 0
        seqs = [
          "AGGGTCACGTAATGCTGATC",
          "GGGTCACGTAATGCTGATCC",
          "GGTCACGTAATGCTGATCCA",
          "GTCACGTAATGCTGATCCAG",
          "TCACGTAATGCTGATCCAGT"
        ]

        instance = GenomeTiler.new
        instance.each_window_in_data(data, 20) do |definition, seq|
          expect(seq).to eq seqs[i]
          i+=1
          break if i == 4
        end
      end

      let(:multiple_data) {
        File.join(
          File.dirname(__FILE__), '..', 'fixtures', 'two_sequence_fragments.fasta'
        )
      }

      context "given the wildcard chromosome filter" do

        it "yields all chromosomes and position-specific subsequence" do
          i = 0
          seqs = [
            "AGGGTCACGTAATGCTGATC",
            "GGGTCACGTAATGCTGATCC",
            "GGTCACGTAATGCTGATCCA",
            "GTCACGTAATGCTGATCCAG",
            "TCACGTAATGCTGATCCAGT",
            "CGGGTCACGTAATGCTGATC",
            "GGGTCACGTAATGCTGATCC",
            "GGTCACGTAATGCTGATCCT",
            "GTCACGTAATGCTGATCCTG",
            "TCACGTAATGCTGATCCTGT"
          ]

          instance = GenomeTiler.new
          results = []
          instance.each_window_in_data(multiple_data, 20, {:filter_chr => ["YHet", "XHet"]}) do |definition, seq|
            results << seq
          end

          expect(results.sort).to eq seqs.sort
        end
      end

      context "given the wildcard chromosome filter" do
        it "yields all chromosomes and position-specific subsequence" do
          i = 0
          seqs = [
            "AGGGTCACGTAATGCTGATC",
            "GGGTCACGTAATGCTGATCC",
            "GGTCACGTAATGCTGATCCA",
            "GTCACGTAATGCTGATCCAG",
            "TCACGTAATGCTGATCCAGT",
            "CGGGTCACGTAATGCTGATC",
            "GGGTCACGTAATGCTGATCC",
            "GGTCACGTAATGCTGATCCT",
            "GTCACGTAATGCTGATCCTG",
            "TCACGTAATGCTGATCCTGT"
          ]

          instance = GenomeTiler.new
          results = []
          instance.each_window_in_data(multiple_data, 20, {:filter_chr => ["*"]}) do |definition, seq|
            results << seq
          end

          expect(results.sort).to eq seqs.sort
        end
      end

      context "given a specific chromosome filter" do
        it "yields a chromosome and position-specific subsequence" do
          i = 0
          seqs = [
            "CGGGTCACGTAATGCTGATC",
            "GGGTCACGTAATGCTGATCC",
            "GGTCACGTAATGCTGATCCT",
            "GTCACGTAATGCTGATCCTG",
            "TCACGTAATGCTGATCCTGT"
          ]

          instance = GenomeTiler.new
          results = []
          instance.each_window_in_data(multiple_data, 20, {:filter_chr => ["YHet"]}) do |definition, seq|
            results << seq
          end
          expect(results.sort).to eq seqs.sort
        end
      end

      it "contains the right chromosome and position in the definition line" do
        i = 0
        definitions = [
          ">yhet_1_20", ">yhet_2_21", ">yhet_3_22", ">yhet_4_23", ">yhet_5_24"
        ]

        instance = GenomeTiler.new
        instance.each_window_in_data(data, 20) do |definition, seq|
          expect(definition).to eq definitions[i]
          i+=1
          break if i == 4
        end
      end

      context "given the option :shifted == 10" do
        it "yields a shifted position-specific subsequence" do
          i = 0
          seqs = [
            "CACGTAATGCTGATCCAGTC", "ACGTAATGCTGATCCAGTCT",
            "CGTAATGCTGATCCAGTCTT", "GTAATGCTGATCCAGTCTTG",
            "TAATGCTGATCCAGTCTTGT"
          ]

          instance = GenomeTiler.new
          instance.each_window_in_data(data, 20, {:shifted => 5}) do |definition, seq|
            expect(seq).to eq seqs[i]
            i+=1
            break if i == 4
          end
        end

        it "contains the correctly shifted chromosome and position in the definition line" do
          i = 0
          definitions = [
            ">yhet_6_25", ">yhet_7_26", ">yhet_8_27", ">yhet_9_28", ">yhet_10_29"
          ]

          instance = GenomeTiler.new
          instance.each_window_in_data(data, 20, {:shifted => 5}) do |definition, seq|
            expect(definition).to eq definitions[i]
            i+=1
            break if i == 4
          end
        end
      end

      describe "#generate_definition" do
        context "given a fasta definition line that contains an ID element" do
          it "generates a definiton string with the ID and start/end position" do
            instance = GenomeTiler.new
            expect(instance.generate_definition(
              {"type"=>"chromosome_arm", "loc"=>"YHet:1..347038", "ID"=>"YHet"},
              5,
              20
            )).to eq ">yhet_6_25"
          end
        end

        describe "#definition_to_fields" do
          context "given a fasta definition line we split it by ';' and map it into a Hash" do
            it "??" do
              instance = GenomeTiler.new
              titles = ["type", "loc", "ID", "dbxref", "MD5", "length", "release", "species"]
              fields = instance.definition_to_fields(">YHet type=chromosome_arm; loc=YHet:1..347038; ID=YHet; dbxref=REFSEQ:NW_001848860,GB:CM000461; MD5=478fbc07ea1b1c03b3d8d04abacf51a5; length=347038; release=r5.51; species=Dmel;")
              expect(titles.sort).to eq fields.keys.sort
            end

            context "given a fasta definition line that misses the ID element" do
              it "throws a GenomeTilerFastaMissingIDError exception" do
                instance = GenomeTiler.new
                expect {
                  instance.definition_to_fields(
                    ">chrU type=chromosome_arm; loc=chrU 1..10049037; MD5=4b72bf19979c8466d5b66acca66f1804; length=10049037; release=r5.55; species=Dmel;"
                  )
                }.to raise_error(
                  GenomeTilerFastaMissingIDError,
                    "Definition line '>chrU type=chromosome_arm; loc=chrU 1..10049037; MD5=4b72bf19979c8466d5b66acca66f1804; length=10049037; release=r5.55; species=Dmel;' misses the ID element."
                )
              end

            end

            context "given an invalid fasta definiton line that contains multiple assignments in one column" do
              it "throws a GenomeTilerFastaMultipleAssignmentsError exception" do
                instance = GenomeTiler.new
                expect {
                  instance.definition_to_fields(
                    ">chrU type=chromosome_arm; loc=chrU 1..10049037; ID=chrU   MD5=4b72bf19979c8466d5b66acca66f1804; length=10049037; release=r5.55; species=Dmel;"
                  )
                }.to raise_error(
                  GenomeTilerFastaMultipleAssignmentsError,
                    "Definition line '>chrU type=chromosome_arm; loc=chrU 1..10049037; ID=chrU   MD5=4b72bf19979c8466d5b66acca66f1804; length=10049037; release=r5.55; species=Dmel;' has multiple assignments in a single column."
                )
              end
            end
          end
        end

      end

    end

    describe "#write_element_to_stream" do
      it "writes the name of an transposon element to a stream" do
        instance = GenomeTiler.new

        output = StringIO.new
        instance.write_element_to_stream(output, ">yhet_1_20", "AGAGAGAGAG")
        output.close
        expect(output.string.split("\n").first).to eq ">yhet_1_20"
      end

      it "writes the sequence of an transposon element to a stream" do
        instance = GenomeTiler.new

        output = StringIO.new
        instance.write_element_to_stream(output, ">yhet_1_20", "AGAGAGAGAG")
        output.close

        expect(output.string.split("\n").last).to eq "AGAGAGAGAG"
      end
    end

  end
end
