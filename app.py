# pipeline_explainer_app/app.py

import streamlit as st
#import openai
from pathlib import Path
import subprocess
import requests
# === CONFIG ===


def get_response(prompt, model="networkjohnny/deepseek-coder-v2-lite-base-q4_k_m-gguf:latest"):
    response = requests.post(
        "http://localhost:11434/api/generate",
        json={"model": model, "prompt": prompt, "stream": False}
    )
    return response.json()["response"]


def explain_code_block(tool, code_block, model):
    prompt = f"Explain the following bioinformatics {tool} code to a junior scientist:\n\n{code_block}. \n\n Just return the explanation, nothing else\n\nExplanation: "

    #response = requests.post(
    #    "https://api-inference.huggingface.co/models/mistralai/Mixtral-8x7B-Instruct-v0.1",
    #    headers={"Authorization": f"Bearer " + st.secrets["HF_API_KEY"]},
    #    json={
    #        "inputs": prompt,
    #        "parameters": {
    #            "temperature": 0.3,
    #            "max_new_tokens": 512,
    #            "return_full_text": False
    #        }
    #    }
    #)
    response = get_response(prompt, model)
    #print(response)

    #result = response.json()
    return response #result[0]["generated_text"].strip() if isinstance(result, list) else "No response or model loading..."

# === HELPER FUNCTIONS ===
def split_script_into_blocks(script_text):
    """
    Split script into logical blocks.
    For Snakemake: split by 'rule'
    For bash: split by function or comment blocks
    """
    blocks = []
    lines = script_text.splitlines()
    current_block = []
    tool = "snakemake" if "rule" in script_text else "bash"
    for line in lines:
        if line.strip().startswith("rule ") or line.strip().startswith("function ") or line.strip().startswith("#"):
            if current_block:
                blocks.append("\n".join(current_block))
                current_block = []
        current_block.append(line)
    if current_block:
        blocks.append("\n".join(current_block))
    return (tool, blocks)

def get_model_names():
    result = subprocess.run(["ollama", "list"], capture_output=True, text=True)
    lines = result.stdout.strip().split('\n')[1:]  # Skip header
    return [line.split()[0] for line in lines]

# === STREAMLIT UI ===
st.set_page_config(page_title="Bioinformatics Pipeline Explainer")
st.title("🧬 Bioinformatics Pipeline Explainer")

uploaded_file = st.file_uploader("Upload a pipeline script (.sh, .smk, .nf)", type=["sh", "smk", "nf", "txt"])
#available_models = 
model = st.selectbox("Choose Model", get_model_names())

if uploaded_file is not None:
    script_text = uploaded_file.read().decode("utf-8")
    st.subheader("📄 Uploaded Script")
    st.code(script_text, language="bash")

    tool, blocks = split_script_into_blocks(script_text)
    st.subheader("🧠 Explanation")

    for i, block in enumerate(blocks):
        with st.expander(f"Step {i+1}"):
            st.code(block, language="bash")
            explanation = explain_code_block(tool, block, model)
            st.markdown(explanation)


