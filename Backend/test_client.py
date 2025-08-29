"""
Test client for the Molecular Simulation FastAPI backend
"""
import requests
import json
import os

# API base URL
BASE_URL = "http://localhost:8001"

def test_upload_and_process():
    """Test file upload and processing"""
    # Check if our test MOL2 file exists
    test_file_path = "117_ideal.mol2"
    
    if not os.path.exists(test_file_path):
        print(f"Test file {test_file_path} not found")
        return None
    
    # Prepare the file for upload
    files = [
        ('files', ('117_ideal.mol2', open(test_file_path, 'rb'), 'chemical/x-mol2'))
    ]
    
    # Test data
    data = {
        'approach': 'both'  # Can be 'classical', 'quantum', or 'both'
    }
    
    try:
        print("Testing file upload and processing...")
        response = requests.post(f"{BASE_URL}/upload-and-process", files=files, data=data)
        
        if response.status_code == 200:
            result = response.json()
            print("✅ Upload and processing successful!")
            print(f"Session ID: {result['session_id']}")
            return result['session_id']
        else:
            print(f"❌ Error: {response.status_code}")
            print(response.text)
            return None
            
    except Exception as e:
        print(f"❌ Connection error: {str(e)}")
        return None
    finally:
        # Close file
        for _, (_, file_obj, _) in files:
            file_obj.close()

def test_get_results(session_id):
    """Test getting results for a session"""
    try:
        print(f"Getting results for session {session_id}...")
        response = requests.get(f"{BASE_URL}/results/{session_id}")
        
        if response.status_code == 200:
            result = response.json()
            print("✅ Results retrieved successfully!")
            print(json.dumps(result, indent=2))
        else:
            print(f"❌ Error getting results: {response.status_code}")
            print(response.text)
            
    except Exception as e:
        print(f"❌ Connection error: {str(e)}")

def test_get_summary(session_id):
    """Test getting summary plots for a session"""
    try:
        print(f"Getting summary for session {session_id}...")
        response = requests.get(f"{BASE_URL}/summary/{session_id}")
        
        if response.status_code == 200:
            result = response.json()
            print("✅ Summary retrieved successfully!")
            
            # Save plot if available
            if 'summary_plots' in result and 'summary_plot' in result['summary_plots']:
                plot_data = result['summary_plots']['summary_plot']
                with open(f"summary_plot_{session_id}.png", "wb") as f:
                    import base64
                    f.write(base64.b64decode(plot_data))
                print(f"📊 Summary plot saved as summary_plot_{session_id}.png")
            
            print(f"Statistics: {json.dumps(result['summary_plots'].get('statistics', {}), indent=2)}")
        else:
            print(f"❌ Error getting summary: {response.status_code}")
            print(response.text)
            
    except Exception as e:
        print(f"❌ Connection error: {str(e)}")

def test_health_check():
    """Test the health check endpoint"""
    try:
        print("Testing health check...")
        response = requests.get(f"{BASE_URL}/health")
        
        if response.status_code == 200:
            result = response.json()
            print("✅ Health check successful!")
            print(f"Status: {result['status']}")
        else:
            print(f"❌ Health check failed: {response.status_code}")
            
    except Exception as e:
        print(f"❌ Connection error: {str(e)}")

if __name__ == "__main__":
    print("🧪 Testing Molecular Simulation API")
    print("="*50)
    
    # Test health check
    test_health_check()
    print()
    
    # Test file upload and processing
    session_id = test_upload_and_process()
    print()
    
    if session_id:
        # Test getting results
        test_get_results(session_id)
        print()
        
        # Test getting summary
        test_get_summary(session_id)
        print()
        
        print(f"🧹 Don't forget to cleanup session: DELETE {BASE_URL}/cleanup/{session_id}")
    
    print("🏁 Testing complete!")
