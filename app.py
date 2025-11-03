import os
import subprocess
from flask import Flask, request, redirect, url_for, render_template, flash, send_file
from werkzeug.utils import secure_filename

# Flask 应用配置
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['OUTPUT_FOLDER'] = 'outputs'
app.config['SECRET_KEY'] = 'your_secret_key'  # 用于安全性

# 确保必要目录存在
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['OUTPUT_FOLDER'], exist_ok=True)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # 处理 BLAST 查询
        if 'query_file' in request.files and 'blast_type' in request.form:
            query_file = request.files['query_file']
            blast_type = request.form['blast_type']
            evalue = request.form['evalue']

            if query_file.filename == '':
                flash('未选择文件')
                return redirect(request.url)
            try:
                float(evalue)  # 确保 evalue 是一个有效的浮点数
            except ValueError:
                flash('无效的 E-value，请输入有效数字')
                return redirect(request.url)

            filename = secure_filename(query_file.filename)
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            query_file.save(file_path)

            # 根据 BLAST 类型选择数据库和命令
            try:
                if blast_type == "nucleotide":
                    db = "lily_db"  # 核酸数据库
                    blast_command = f"blastn -query {file_path} -db {db} -out {file_path}.raw -evalue {evalue} -outfmt 6 -num_threads 10"
                    output_file = os.path.join(app.config['OUTPUT_FOLDER'], f"{filename}.results")
                    conversion_command = f"python3 conversion_blast_results.py {file_path}.raw {output_file}"  # 核酸结果转换
                elif blast_type == "protein":
                    db = "lily_prot_db"  # 蛋白质数据库
                    output_file=os.path.join(app.config['OUTPUT_FOLDER'],f"{filename}.prot.results")
                    blast_command = f"blastp -query {file_path} -db {db} -out {output_file} -evalue {evalue} -outfmt 6 -num_threads 10"
                    #output_file = os.path.join(app.config['OUTPUT_FOLDER'], f"{filename}.prot.results")
                    conversion_command = None  # 蛋白查询无需结果转换
                else:
                    flash("无效的 BLAST 类型")
                    return redirect(request.url)

                print(f"Running BLAST command: {blast_command}")
                subprocess.run(blast_command, shell=True, check=True, timeout=300)

                # 仅在核酸 BLAST 查询时运行结果转换脚本
                if blast_type == "nucleotide" and conversion_command:
                    print(f"Running conversion command: {conversion_command}")
                    subprocess.run(conversion_command, shell=True, check=True)

                # 提供下载链接
                return redirect(url_for('results', result_file=os.path.basename(output_file)))
            except subprocess.CalledProcessError as e:
                flash(f"BLAST 或脚本执行失败: {e}")
                return redirect(request.url)
            except subprocess.TimeoutExpired:
                flash("BLAST 执行超时，请稍后重试")
                return redirect(request.url)

        # 处理染色体范围查询 (功能2)
        elif 'chromosome' in request.form and 'start' in request.form and 'end' in request.form:
            chromosome = request.form['chromosome']
            start = request.form['start']
            end = request.form['end']

            # 确保起点和终点是整数
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                flash("起点和终点必须是整数")
                return redirect(request.url)

            # 生成临时 BED 文件
            bed_file = os.path.join(app.config['UPLOAD_FOLDER'], 'temp.bed')
            with open(bed_file, 'w') as bed:
                bed.write(f"{chromosome}\t{start-1}\t{end}\n")  # BED 文件起点为 0 基

            # 定义输出文件路径
            genome_file = "Ldavi.chr.fasta"  # 假设这是你的基因组文件
            output_file = os.path.join(app.config['OUTPUT_FOLDER'], f"{chromosome}_{start}_{end}.fa")

            # 运行 seqkit subseq 命令
            seqkit_command = f"seqkit subseq --bed {bed_file} {genome_file} -o {output_file}"
            try:
                print(f"Running seqkit command: {seqkit_command}")
                subprocess.run(seqkit_command, shell=True, check=True)

                # 删除临时 BED 文件
                os.remove(bed_file)

                # 提供下载链接
                return redirect(url_for('results', result_file=os.path.basename(output_file)))
            except subprocess.CalledProcessError as e:
                flash(f"序列提取失败: {e}")
                return redirect(request.url)

    return render_template('index.html')


@app.route('/results/<result_file>')
def results(result_file):
    result_path = os.path.join(app.config['OUTPUT_FOLDER'], result_file)
    if not os.path.exists(result_path):
        return "结果文件不存在", 404
    return render_template('results.html', result_file=result_file)


@app.route('/download/<result_file>')
def download(result_file):
    result_path = os.path.join(app.config['OUTPUT_FOLDER'], result_file)
    if not os.path.exists(result_path):
        return "文件不存在", 404
    return send_file(result_path, as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
