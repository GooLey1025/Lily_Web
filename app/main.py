from fastapi import FastAPI, Request, UploadFile, File, Form
from fastapi.responses import HTMLResponse, JSONResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
import os
import uuid
import shutil
import subprocess
from typing import Dict, Any


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(BASE_DIR)
DATA_DIR = os.path.join(ROOT_DIR, "data")
TEMPLATES_DIR = os.path.join(ROOT_DIR, "app", "templates")
STATIC_DIR = os.path.join(ROOT_DIR, "app", "static")

os.makedirs(DATA_DIR, exist_ok=True)

app = FastAPI()
templates = Jinja2Templates(directory=TEMPLATES_DIR)
app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


def _save_upload(file: UploadFile, dest_path: str) -> None:
    with open(dest_path, "wb") as out:
        shutil.copyfileobj(file.file, out)


def _run(cmd: list[str], cwd: str | None = None, timeout: int | None = 0) -> None:
    subprocess.run(cmd, cwd=cwd, check=True)


def _parse_outfmt6_to_json(lines: list[str]) -> list[Dict[str, Any]]:
    results = []
    for ln in lines:
        ln = ln.strip()
        if not ln:
            continue
        cols = ln.split("\t")
        if len(cols) < 12:
            # allow shorter if needed for robustness
            pass
        # Standard outfmt 6 columns
        entry = {
            "qseqid": cols[0],
            "sseqid": cols[1],
            "pident": float(cols[2]),
            "length": int(cols[3]),
            "mismatch": int(cols[4]),
            "gapopen": int(cols[5]),
            "qstart": int(cols[6]),
            "qend": int(cols[7]),
            "sstart": int(cols[8]),
            "send": int(cols[9]),
            "evalue": float(cols[10]) if len(cols) > 10 else None,
            "bitscore": float(cols[11]) if len(cols) > 11 else None,
            "raw": ln,
        }
        results.append(entry)
    return results


@app.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})


@app.post("/blastn")
async def run_blastn(
    file: UploadFile = File(...),
):
    try:
        job_id = str(uuid.uuid4())
        job_dir = os.path.join(DATA_DIR, job_id)
        os.makedirs(job_dir, exist_ok=True)

        query_path = os.path.join(job_dir, "query.fasta")
        _save_upload(file, query_path)

        raw_path = os.path.join(job_dir, "blast.raw.txt")
        result_path = os.path.join(job_dir, "blast.results.txt")

        _run([
            "blastn", "-query", query_path, "-db", "lily_db",
            "-out", raw_path, "-evalue", "1e-5", "-outfmt", "6", "-num_threads", "10",
        ])

        # Convert to absolute coordinates; only expose converted file
        converter = os.path.join(ROOT_DIR, "conversion_blast_results.py")
        _run(["python3", converter, raw_path, result_path])

        with open(result_path, "r", encoding="utf-8") as f:
            json_rows = _parse_outfmt6_to_json(f.readlines())

        return JSONResponse({
            "job_id": job_id,
            "rows": json_rows,
            "download": {
                # raw not exposed by requirement
                "results": f"/download/{job_id}/results",
            },
        })
    except subprocess.CalledProcessError as e:
        return JSONResponse(
            {"error": f"BLASTN execution failed: {str(e)}"},
            status_code=500
        )
    except Exception as e:
        return JSONResponse(
            {"error": f"BLASTN error: {str(e)}"},
            status_code=500
        )


@app.post("/blastp")
async def run_blastp(
    file: UploadFile = File(...),
):
    try:
        job_id = str(uuid.uuid4())
        job_dir = os.path.join(DATA_DIR, job_id)
        os.makedirs(job_dir, exist_ok=True)

        query_path = os.path.join(job_dir, "query.fasta")
        _save_upload(file, query_path)

        raw_path = os.path.join(job_dir, "blastp.raw.txt")
        result_path = os.path.join(job_dir, "blastp.results.txt")

        _run([
            "blastp", "-query", query_path, "-db", "lily_prot_db",
            "-out", raw_path, "-evalue", "1e-5", "-outfmt", "6", "-num_threads", "10",
        ])

        converter = os.path.join(ROOT_DIR, "conversion_blast_results.py")
        _run(["python3", converter, raw_path, result_path])

        with open(result_path, "r", encoding="utf-8") as f:
            json_rows = _parse_outfmt6_to_json(f.readlines())

        return JSONResponse({
            "job_id": job_id,
            "rows": json_rows,
            "download": {
                "results": f"/download/{job_id}/results",
            },
        })
    except subprocess.CalledProcessError as e:
        return JSONResponse(
            {"error": f"BLASTP execution failed: {str(e)}"},
            status_code=500
        )
    except Exception as e:
        return JSONResponse(
            {"error": f"BLASTP error: {str(e)}"},
            status_code=500
        )


@app.post("/subseq")
async def run_subseq(
    chromosome: str = Form(...),
    start: int = Form(...),
    end: int = Form(...),
):
    try:
        job_id = str(uuid.uuid4())
        job_dir = os.path.join(DATA_DIR, job_id)
        os.makedirs(job_dir, exist_ok=True)

        # Generate BED file from 1-based coordinates (convert to 0-based)
        bed_path = os.path.join(job_dir, "user.bed")
        bed_start = start - 1  # Convert to 0-based
        bed_end = end  # End remains the same for half-open interval
        
        with open(bed_path, "w") as bed_file:
            bed_file.write(f"{chromosome}\t{bed_start}\t{bed_end}\n")

        out_fa = os.path.join(job_dir, "user.fa")

        # seqkit writes to stdout; capture to file
        with open(out_fa, "wb") as out:
            subprocess.run([
                "seqkit", "subseq", "--bed", bed_path, "Ldavi.chr.fasta"
            ], check=True, stdout=out, stderr=subprocess.PIPE)

        return JSONResponse({
            "job_id": job_id,
            "download": {
                "fasta": f"/download/{job_id}/subseq",
            },
        })
    except subprocess.CalledProcessError as e:
        error_msg = e.stderr.decode('utf-8') if e.stderr else str(e)
        return JSONResponse(
            {"error": f"Sequence extraction failed: {error_msg}"},
            status_code=500
        )
    except Exception as e:
        return JSONResponse(
            {"error": f"Subseq error: {str(e)}"},
            status_code=500
        )


@app.get("/result/{job_id}")
async def get_result(job_id: str):
    job_dir = os.path.join(DATA_DIR, job_id)
    if not os.path.isdir(job_dir):
        return JSONResponse({"error": "job not found"}, status_code=404)

    # Prefer blast results if present
    for fname in ["blast.results.txt", "blastp.results.txt"]:
        fpath = os.path.join(job_dir, fname)
        if os.path.exists(fpath):
            with open(fpath, "r", encoding="utf-8") as f:
                return JSONResponse({"job_id": job_id, "rows": _parse_outfmt6_to_json(f.readlines())})

    return JSONResponse({"job_id": job_id, "rows": []})


@app.get("/download/{job_id}/{which}")
async def download(job_id: str, which: str):
    job_dir = os.path.join(DATA_DIR, job_id)
    if which == "results":
        # Prefer nucleotide, then protein
        for fname in ["blast.results.txt", "blastp.results.txt"]:
            fpath = os.path.join(job_dir, fname)
            if os.path.exists(fpath):
                return FileResponse(fpath, filename=fname, media_type="text/plain")
        return JSONResponse({"error": "results not found"}, status_code=404)
    if which == "subseq":
        fpath = os.path.join(job_dir, "user.fa")
        if os.path.exists(fpath):
            return FileResponse(fpath, filename="user.fa", media_type="application/octet-stream")
        return JSONResponse({"error": "fasta not found"}, status_code=404)
    return JSONResponse({"error": "unsupported"}, status_code=400)


