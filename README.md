# Efficient DNA Sequence Compression: A Comparative Study of Huffman and Arithmetic Coding with Data Structures

🔬 **Project**: Comparative Study of Huffman and Arithmetic Coding Techniques for Efficient DNA Sequence Compression  
📄 **Summary**: Designed and implemented a modular system to compress DNA sequences using two techniques: **Huffman Encoding** (📊 prefix-based encoding) and **Arithmetic Coding** (🎲 probability-based encoding).

---

## ⚙️ Key Features

- **Huffman Coding**:  
   - Constructed a **Binary Search Tree (BST)** for efficient prefix-based encoding.  
   - Suitable for sequences with non-uniform symbol frequencies, ensuring optimal performance.  
   
- **Arithmetic Coding**:  
   - Utilized **probability intervals** for finer encoding granularity.  
   - Offers better compression for sequences with highly repetitive or predictable patterns.

- **Memory Optimization**:  
   - Both techniques demonstrated **significant memory savings** when applied to DNA sequences, enhancing data storage and transmission.

---

## 📈 Results

- **Huffman Coding**:  
   - **Compression Ratio**: 0.33 (significantly reduced data size).  
   - **Processing Speed**: Faster encoding/decoding due to simple tree structure.

- **Arithmetic Coding**:  
   - **Compression Ratio**: 0.051 (achieved higher compression granularity).  
   - Provides superior efficiency when dealing with highly repetitive data, like DNA sequences.

---

## 📚 Impact

- **Enhanced Data Storage**: The project addresses challenges in **genomic data storage** by optimizing the size of DNA sequence files.
- **Improved Data Transmission**: More efficient compression leads to **faster data transfer** for bioinformatics applications, reducing storage and computational costs.
- **Bioinformatics**: The compression methods explored can potentially support larger-scale genomic analysis and data sharing, making them valuable for the field.

---

## 💡 Tech Stack

- **Programming Language**: Python  
- **Key Concepts**: Data Structures (BST, Probability Intervals), Algorithm Design, Compression Techniques  
- **Compression Algorithms**: Huffman Encoding, Arithmetic Coding  
