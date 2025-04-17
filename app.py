# pipeline_explainer_app/app.py

import streamlit as st
import openai
from pathlib import Path

# === CONFIG ===
import requests

def explain_code_block(code_block):
    prompt = f"Explain the following bioinformatics code to a junior scientist:\n\n{code_block}\n\nExplanation:"

    response = requests.post(
        "https://api-inference.huggingface.co/models/mistralai/Mixtral-8x7B-Instruct-v0.1",
        headers={"Authorization": f"Bearer " + st.secrets["HF_API_KEY"]},
        json={
            "inputs": prompt,
            "parameters": {
                "temperature": 0.3,
                "max_new_tokens": 512,
                "return_full_text": False
            }
        }
    )
    print(response.json())

    result = response.json()
    return result[0]["generated_text"].strip() if isinstance(result, list) else "No response or model loading..."

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
    for line in lines:
        if line.strip().startswith("rule ") or line.strip().startswith("function ") or line.strip().startswith("#"):
            if current_block:
                blocks.append("\n".join(current_block))
                current_block = []
        current_block.append(line)
    if current_block:
        blocks.append("\n".join(current_block))
    return blocks


# === STREAMLIT UI ===
st.set_page_config(page_title="Bioinformatics Pipeline Explainer")
st.title("🧬 Bioinformatics Pipeline Explainer")

uploaded_file = st.file_uploader("Upload a pipeline script (.sh, .smk, .nf)", type=["sh", "smk", "nf", "txt"])

if uploaded_file is not None:
    script_text = uploaded_file.read().decode("utf-8")
    st.subheader("📄 Uploaded Script")
    st.code(script_text, language="bash")

    blocks = split_script_into_blocks(script_text)
    st.subheader("🧠 Explanation")

    for i, block in enumerate(blocks):
        with st.expander(f"Step {i+1}"):
            st.code(block, language="bash")
            explanation = explain_code_block(block)
            st.markdown(explanation)

