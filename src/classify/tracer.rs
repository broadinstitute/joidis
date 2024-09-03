use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use log::warn;
use crate::data::Meta;
use crate::error::Error;
use crate::sample::sampler::Tracer;
use crate::sample::var_stats::VarStats;
use crate::sample::vars::Vars;
use crate::util::joiner::{Joiner, write_iter_io};


struct Writer {
    buffer: Vec<f64>,
    writer: Result<BufWriter<File>, Error>
}

pub(crate) struct ClassifyTracer {
    e_writers: Vec<Writer>,
    t_writers: Vec<Writer>,
    conv_writer: Result<BufWriter<File>, Error>,
    n_effective_writer: Result<BufWriter<File>, Error>,
}

impl Writer {
    fn new(file_name: &str, n_chains: usize) -> Writer {
        let buffer = vec![f64::NAN; n_chains];
        let writer = try_writer(file_name);
        Writer { buffer, writer }
    }
    fn set_value(&mut self, i_chain: usize, value: f64) {
        self.buffer[i_chain] = value;
    }
    fn write_values(&mut self) {
        if let Ok(ref mut writer) = self.writer {
            if let Err(error) = writeln!(writer, "{}", Joiner::new("\t", &self.buffer)) {
                warn!("Could not write trace: {}", error)
            }
        }
    }
}

impl ClassifyTracer {
    pub(crate) fn new(meta: &Meta, out_file_name: &str, var_id: &str, n_endos: usize,
                      n_chains: usize) -> ClassifyTracer {
        let e_writers = (0..n_endos).map(|i_endo| {
            let var_name = format!("E_{}", i_endo);
            let file_name = format!("{}_{}_trace_{}", out_file_name, var_id, var_name);
            Writer::new(&file_name, n_chains)
        }).collect::<Vec<_>>();
        let t_writers = (0..meta.n_traits()).map(|i_trait| {
            let var_name = format!("T_{}", i_trait);
            let file_name = format!("{}_{}_trace_{}", out_file_name, var_id, var_name);
            Writer::new(&file_name, n_chains)
        }).collect::<Vec<_>>();
        let conv_file_name = format!("{}_{}_trace_conv", out_file_name, var_id);
        let mut conv_writer = try_writer(&conv_file_name);
        if let Ok(conv_writer) = &mut conv_writer {
            let variable_names = Vars::variable_names(meta);
            warn_if_err(write_iter_io(conv_writer, variable_names, "\t"));
            warn_if_err(writeln!(conv_writer))
        }
        let n_effective_file_name = format!("{}_{}_trace_n_eff", out_file_name, var_id);
        let n_effective_writer = try_writer(&n_effective_file_name);
        ClassifyTracer { e_writers, t_writers, conv_writer, n_effective_writer }
    }
}

fn warn_if_err<T>(result: std::io::Result<T>) {
    if let Err(error) = result {
        warn!("{}", error)
    }
}

fn try_writer(file_name: &str) -> Result<BufWriter<File>, Error> {
    match File::create(file_name) {
        Ok(file) => { Ok(BufWriter::new(file)) }
        Err(error) => { Err(Error::from(error)) }
    }
}

impl Tracer for ClassifyTracer {
    fn trace_e(&mut self, i_endo: usize, e: f64, i_chain: usize) {
        self.e_writers[i_endo].set_value(i_chain, e);
    }

    fn trace_t(&mut self, i_trait: usize, t: f64, i_chain: usize) {
        self.t_writers[i_trait].set_value(i_chain, t);
    }
    fn write_values(&mut self) {
        for writer in self.e_writers.iter_mut() {
            writer.write_values()
        }
        for writer in self.t_writers.iter_mut() {
            writer.write_values()
        }
    }

    fn trace_convergence(&mut self, var_stats_list: &[VarStats]) {
        if let Ok(ref mut writer) = self.conv_writer {
            let convergences =
                VarStats::calculate_convergences(var_stats_list);
            warn_if_err(write_iter_io(writer, convergences, "\t"));
            warn_if_err(writeln!(writer));
        }
    }

    fn trace_n_effective(&mut self, var_stats_list: &[VarStats]) {
        if let Ok(ref mut writer) = self.n_effective_writer {
            let n_effectives =
                var_stats_list.iter().map(|var_stats: &VarStats|
                    var_stats.effective_sample_size()
                );
            warn_if_err(write_iter_io(writer, n_effectives, "\t"));
            warn_if_err(writeln!(writer));
        }
    }
}
