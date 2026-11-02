# Documentation helper only. It does not fit a model or choose a scientific prior.
require 'json'
abort 'Usage: ruby prepare_config.rb MODEL_JSON OUTPUT_JSON [PRIOR_TEMPLATE_JSON]' unless (2..3).cover?(ARGV.length)
model_path, output_path, template_path = ARGV
abort 'Refusing existing config output' if File.exist?(output_path) || File.symlink?(output_path)
template_path ||= 'fixtures/synthetic/hybrid.json'
model = JSON.parse(File.read(model_path))
abort 'A successful qualified legacy model is required' unless model['mode'] == 'legacy-compatible' && model.dig('fit', 'success') == true
parameters = model.fetch('fit').fetch('evaluation').fetch('parameters')
config = JSON.parse(File.read(template_path))
config.fetch('noise')['parameters'] = [parameters.fetch(0), 10.0 ** parameters.fetch(1)]
config['noise']['legacy_model_file'] = File.absolute_path(model_path)
config['noise']['source'] = 'qualified_legacy_fit'
config['noise']['source_reference'] = 'Successful Core fit; response prior supplied by template'
File.open(output_path, File::WRONLY | File::CREAT | File::EXCL) { |f| f.write(JSON.pretty_generate(config) + "\n") }
puts "Created #{output_path}; review response_prior and effect_delta before scientific use"
