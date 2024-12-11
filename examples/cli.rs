use clap::Parser;
use lagrit_bindings::{InitMode, LaGriT, LagritError};
use lagrit_sys::fc_control_command_lg;

#[derive(Parser)]
struct Opts {
    #[arg(short, long, default_value = "lagrit.log")]
    log: String,
    #[arg(short, long, default_value = "lagrit.out")]
    out: String,
    #[arg(value_enum, short, long, default_value = "noisy")]
    mode: InitMode,
}

fn main() -> Result<(), LagritError> {
    let cli_opts: Opts = Opts::parse();

    LaGriT::new(
        cli_opts.mode,
        Some(&cli_opts.log),
        Some(&cli_opts.out),
        None,
    )?;

    let mut status: i32 = 0;
    unsafe { fc_control_command_lg(&mut status) };

    Ok(())
}
