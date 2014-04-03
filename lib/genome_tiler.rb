require 'stringio'
require 'bio'

class GenomeTiler
  def split_into_windows(data, output, window_size, options)
    items = 0

    each_sequence_in_data(data, window_size) do |definition, sequence|
      write_element_to_stream(output, definition, sequence)
      items += 1
    end

    items
  end

  def each_sequence_in_data(data, window_size)
    Bio::FlatFile.auto(data) do |ff|
      ff.each do |entry|
        definition = entry.definition
        sequence = entry.seq
        (0..sequence.length-window_size).each do |i|
          yield generate_definition(definition, i, window_size), sequence[i..(i+window_size-1)]
        end
      end
    end
  end

  def generate_definition(definition, position, window_size)
    fields = definition_to_fields(definition)
    ">#{fields["ID"].downcase}_#{position}_#{position+window_size}"
  end

  def definition_to_fields(definition)
    definition.split(";").map { |s| s.strip.split("=") }.reduce({}) { |o, n| o.merge!(n[0] => n[1]) }
  end

  def write_element_to_stream(output, definition, sequence)
    output.puts("#{definition}\n#{sequence}\n")
  end
end

