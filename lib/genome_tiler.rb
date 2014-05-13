require 'command'
require 'stringio'
require 'bio'

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
        definition = entry.definition
        sequence = entry.seq
        start_pos = options.fetch(:shifted) { 0 }

        (start_pos..sequence.length-window_size).each do |i|
          yield generate_definition(definition, i, window_size), sequence[i..(i+window_size-1)]
        end
      end
    end
  end

  def generate_definition(definition, position, window_size)
    fields = definition_to_fields(definition)
    raise "Definition line #{definition} misses ID element [#{window_size}]" unless fields.has_key?("ID")
    ">#{fields["ID"].downcase}_#{position + 1}_#{position+window_size}"
  end

  def definition_to_fields(definition)
    definition.split(";").map { |s| s.strip.split("=") }.reduce({}) { |o, n| o.merge!(n[0] => n[1]) }
  end

  def write_element_to_stream(output, definition, sequence)
    output.puts("#{definition}\n#{sequence}\n")
  end
end

