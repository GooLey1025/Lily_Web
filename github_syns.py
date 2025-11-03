import requests
import time
import base64

# GitHub 配置信息
GITHUB_TOKEN = ""  # 你的 GitHub Token
REPO_OWNER = "GooLey1025"  # 仓库拥有者
REPO_NAME = "Lily_web_ip_syns"  # 仓库名称
FILE_PATH = "README.md"  # 仓库中的文件路径

def get_ngrok_url():
    """获取当前 ngrok 的外网地址"""
    try:
        response = requests.get("http://127.0.0.1:4040/api/tunnels")
        response.raise_for_status()
        tunnels = response.json().get("tunnels", [])
        for tunnel in tunnels:
            if tunnel.get("proto") == "https":
                return tunnel.get("public_url")
    except Exception as e:
        print(f"Error fetching ngrok URL: {e}")
    return None

def upload_to_github(url):
    """将 ngrok URL 上传到 GitHub 文件"""
    # 获取文件的当前内容和 SHA 值
    api_url = f"https://api.github.com/repos/{REPO_OWNER}/{REPO_NAME}/contents/{FILE_PATH}"
    headers = {
        "Authorization": f"token {GITHUB_TOKEN}",
        "Accept": "application/vnd.github.v3+json",
    }

    try:
        # 检查文件是否已经存在
        response = requests.get(api_url, headers=headers)
        if response.status_code == 200:
            file_data = response.json()
            sha = file_data["sha"]
        else:
            sha = None  # 文件不存在

        # 准备要上传的内容
		# 准备要上传的内容
        content=(
            f"Lily_web URL: {url}\n\n"
            f"Updated at: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n"
            "If you cannot open this URL, the server might be down.\n\n"
            "Please contact: \n\n"
            "goley04@foxmail.com Or Wechat:GL_yeaH"
        )
        encoded_content = base64.b64encode(content.encode("utf-8")).decode("utf-8")

        # 构建上传数据
        data = {
            "message": "Update ngrok URL",
            "content": encoded_content,
        }
        if sha:
            data["sha"] = sha

        # 上传文件到 GitHub
        response = requests.put(api_url, json=data, headers=headers)
        response.raise_for_status()
        print("Ngrok URL successfully uploaded to GitHub.")
    except Exception as e:
        print(f"Error uploading to GitHub: {e}")

if __name__ == "__main__":
    time.sleep(10)  # 等待 ngrok 启动
    url = get_ngrok_url()
    if url:
        upload_to_github(url)
    else:
        print("Failed to retrieve ngrok URL.")
