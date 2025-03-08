# QPharma
Quantum AI Exploration for Drug Discovery

# Quantum AI Exploration for Drug Discovery

This project explores **Quantum Computing** for drug discovery using **Amazon Braket** and integrates it with a **React.js frontend** for an intuitive user experience. The system allows users to upload **.mol2** files, process them via AWS, and retrieve results for analysis.

---

## 🚀 Features

- **Upload .mol2 Files**: Users can upload molecular structure files.
- **Quantum Processing**: AWS Lambda & Amazon Braket handle computations.
- **Results Display**: Fetch and display analyzed JSON results.
- **Secure & Scalable**: Utilizes AWS services (S3, API Gateway, Lambda).
- **User-Friendly UI**: Built with React.js & TailwindCSS.

---

## 🛠️ Tech Stack



### **Backend**
- **Node.js & Express**
- **AWS Lambda (Python)**
- **AWS API Gateway**
- **AWS S3 (Storage)**
- **Amazon Braket (Quantum Processing)**

---

## 📂 Project Structure

Quantum-Drug-Discovery/
│── frontend/               # React Frontend
│   ├── src/
│   │   ├── components/     # UI Components
│   │   ├── App.jsx         # Main React Component
│   │   ├── index.css       # Styles
│   ├── package.json        # Dependencies
│── backend/                # Node.js Backend API
│   ├── server.js           # Express Server
│   ├── package.json        # Dependencies
│── README.md               # Documentation
