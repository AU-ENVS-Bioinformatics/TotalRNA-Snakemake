import datetime
from pathlib import Path
import subprocess
from typing import Generator
import typer
app = typer.Typer()

def create_temporary_dir(prefix: str)-> str:
    suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
    return f"{prefix}_{suffix}"
def echo_message(message: str) -> None:
    typer.secho(message)
def echo_success(message: str) -> None:
    typer.secho(message, fg=typer.colors.GREEN)
    typer.secho("")
def echo_error_and_exit(message: str) -> None:
    typer.secho(message, fg=typer.colors.RED)
    raise typer.Exit(code=1)
def debug_command(commmand: list) -> None:
    typer.secho(" ".join(commmand), fg=typer.colors.MAGENTA)
def construct_translation_command(infile: Path, outfile: Path, orf: int) -> "list[str]":
    """
    Construct translation command using transeq.
    """
    return [
        "transeq", "-sformat", "pearson",
        "-clean", "-frame", str(orf),
        "-sequence", str(infile),
         "-outseq", str(outfile)
        ]
def construct_splitting_command(infile: Path, splitsize: int) -> "list[str]":
    """
    Construct splitting command using pyfasta.
    """
    return ["pyfasta", "split", "-n", str(splitsize), str(infile)]

def get_splitted_fasta_files(translated_fasta: Path) -> list[Path]:
    x = [fasta for fasta in translated_fasta.parent.glob('*.fasta')]
    if len(x) > 1:
        x.remove(translated_fasta)
    return x
def construct_split_outfile(fasta: Path, database: Path):
        return str(fasta.parent / f"{fasta.name}_{database.name}_result.tsv")
def construct_batchsword_command(
    infiles: Generator[Path, None, None],
    database: Path,
    threads: int
    ) -> "list[str]":
    """
    Construct batchsword command using sword.
    """
    def inner(fasta: Path, database: Path):
        return [
            "sword", "-i", str(fasta), "-t", str(threads),
            "-o", construct_split_outfile(fasta, database),
            "-f", "bm9",  "-j", str(database), "-c", "30000"
        ]
    return [inner(fasta, database) for fasta in infiles]
def concatenate_files(infiles: list[Path], outfile: Path):
    outfile.parent.mkdir(parents=True, exist_ok=True)
    merging = " ".join(infiles)
    subprocess.run(f"cat {merging} > {str(outfile)}" ,shell= True, check=True)
@app.command()
def align_contigs_to_database(
    infile: Path = typer.Argument(..., help="Fasta file of assembled contigs, output from Trinity"),
    outfile: Path = typer.Argument(..., help="Output tsv file"),
    database: Path = typer.Argument(..., help="Path to database"),
    splitsize: int = typer.Option(1, "-s", "--splitsize", help="Number of parts Fasta file to be splitted in"),
    orf: int = typer.Option(1, "-n", "--ORFs", min=1, max=6,),
    threads: int = typer.Option(1, "-t", "--threads"),
    remove: bool = typer.Option(False, "--remove")
) -> None:
    """
    This script will use SWORD to align the assembled contigs from
    previous step against database of choice from following options
    """
    """
    This script will use SWORD to align the assembled contigs from
    previous step against database of choice from following options
    """
    #Find databases
    # Create temporary directory
    temp_directory = outfile.parent / create_temporary_dir("temp")
    temp_directory.mkdir(exist_ok=True,parents=True)
    echo_success(f"The directory for temporary results was created: {str(temp_directory)}")
    # Translate sequence into aminoacids
    echo_message(f"Translating the query sequence to proteins in {str(orf)} ORF(s)")
    translated_sequences_file = temp_directory / f"translated_{infile.name}"
    translation_command = construct_translation_command(
        infile, translated_sequences_file, orf
        )
    debug_command(translation_command)
    subprocess.run(translation_command, check=True)
    echo_success("Translation done")
    # Split fasta files
    if splitsize > 1:
        echo_message(f"Splitting translated fasta file into {splitsize} to use less memory")
        splitting_command = construct_splitting_command(translated_sequences_file, splitsize)
        debug_command(splitting_command)
        subprocess.run(splitting_command, check=True)
        echo_success("Splitting done")  
    # Align to databases
    splitted_fasta = list(get_splitted_fasta_files(translated_sequences_file))
    batchsword_commands = construct_batchsword_command(splitted_fasta, database, threads)
    for command in batchsword_commands:
        debug_command(command)
        subprocess.run(command, check=True) 
    merging_files = [construct_split_outfile(fasta, database) for fasta in splitted_fasta]
    concatenate_files(merging_files, outfile)
    if remove:
        subprocess.run(["rm", "-rf", str(temp_directory)], check=True)

if __name__ == "__main__":
    app()