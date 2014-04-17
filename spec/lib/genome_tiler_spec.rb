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
        expect { |b| instance.each_window_in_data(data, 20, &b) }.to yield_control.exactly(31).times
      end

      it "yields a position-specific subsequence" do
        i = 0
        seqs = ["AGGGTCACGTAATGCTGATC", "GGGTCACGTAATGCTGATCC", "GGTCACGTAATGCTGATCCA", "GTCACGTAATGCTGATCCAG", "TCACGTAATGCTGATCCAGT"]

        instance = GenomeTiler.new
        instance.each_window_in_data(data, 20) do |definition, seq|
          expect(seq).to eq seqs[i]
          i+=1
          break if i == 4
        end
      end

      it "contains the right chromosome and position in the definition line" do
        i = 0
        definitions = [">yhet_0_20", ">yhet_1_21", ">yhet_2_22", ">yhet_3_23", ">yhet_4_24"]

        instance = GenomeTiler.new
        instance.each_window_in_data(data, 20) do |definition, seq|
          expect(definition).to eq definitions[i]
          i+=1
          break if i == 4
        end
      end

      describe "#generate_definition" do
        context "given a fasta definition line that contains an ID element" do
          it "generates a definiton string with the ID and start/end position" do
            instance = GenomeTiler.new
            expect(instance.generate_definition(">YHet type=chromosome_arm; loc=YHet:1..347038; ID=YHet;", 5, 20)).to eq ">yhet_5_25"
          end
        end

        describe "#definition_to_fields" do
          context "given a fasta definition line we split it by ';' and map it into a Hash" do
            it "??" do
              instance = GenomeTiler.new
              titles = [">YHet type", "loc", "ID", "dbxref", "MD5", "length", "release", "species"]
              fields = instance.definition_to_fields(">YHet type=chromosome_arm; loc=YHet:1..347038; ID=YHet; dbxref=REFSEQ:NW_001848860,GB:CM000461; MD5=478fbc07ea1b1c03b3d8d04abacf51a5; length=347038; release=r5.51; species=Dmel;")
              fields.keys.each do |key|
                expect(titles.include?(key)).to be_true
                titles.delete(key)
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
        instance.write_element_to_stream(output, ">yhet_0_20", "AGAGAGAGAG")
        output.close
        expect(output.string.split("\n").first).to eq ">yhet_0_20"
      end

      it "writes the sequence of an transposon element to a stream" do
        instance = GenomeTiler.new

        output = StringIO.new
        instance.write_element_to_stream(output, ">yhet_0_20", "AGAGAGAGAG")
        output.close

        expect(output.string.split("\n").last).to eq "AGAGAGAGAG"
      end
    end

  end
end
