require 'command'
require 'stringio'
require 'bio'

class GenomeTilerError < StandardError; end;
class GenomeTilerFastaMissingIDError < GenomeTilerError; end;
class GenomeTilerFastaMultipleAssignmentsError < GenomeTilerError; end;

class GenomeTiler
  def split_into_windows(data, output, window_size, options)
    items = 0

    each_window_in_data(data, window_size, options) do |definition, sequence|
      write_element_to_stream(output, definition, sequence)
      items += 1
    end

    items
  end

  def each_sequence(path_to_input)
    Bio::FlatFile.auto(path_to_input) do |ff|
      ff.each do |entry|
        definition = entry.definition
        sequence = entry.seq

        yield definition, sequence
      end
    end
  end

  def reverse_complement(path_to_input, output_stream)
    each_sequence(path_to_input) do |name, sequence|
      s = reverse_complement_sequence(sequence)
      output_stream.puts(">#{name}\n#{s}\n")
    end
  end

  def reverse_complement_sequence(sequence)
    Bio::Sequence::NA.new(sequence).reverse_complement
  end

  def each_window_in_data(data, window_size, options = {})
    Bio::FlatFile.auto(data) do |ff|
      ff.each do |entry|
        definition_fields = definition_to_fields(entry.definition)

        sequence = entry.seq
        start_pos = options.fetch(:shifted) { 0 }

        (start_pos..sequence.length-window_size).each do |i|
          yield generate_definition(definition_fields, i, window_size), sequence[i..(i+window_size-1)]
        end
      end
    end
  end

  def generate_definition(fields, position, window_size)
    ">#{fields["ID"].downcase}_#{position + 1}_#{position+window_size}"
  end

  def definition_to_fields(definition)

    ## first lets get rid of the '>chrom_name' prefix, it should be available
    ## inside the ID column anyway
    definition_string = definition.gsub(/>\w+/, "").strip

    columns = definition_string.split(";")
    fields = columns.map { |s| s.strip.split("=") }

    raise GenomeTilerFastaMultipleAssignmentsError,
      "Definition line '#{definition}' has multiple assignments in a " \
      "single column." unless fields.select { |f| f.count != 2 }.empty?

    fields = fields.reduce({}) { |o, n| o.merge!(n[0] => n[1]) }

    raise GenomeTilerFastaMissingIDError,
      "Definition line '#{definition}' misses the ID element." unless fields.has_key?("ID")

    fields
  end

  def write_element_to_stream(output, definition, sequence)
    output.puts("#{definition}\n#{sequence}\n")
  end
end

