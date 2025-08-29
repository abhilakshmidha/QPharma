import urllib.request
import json

# Simple test without external dependencies
url = "http://localhost:8001/upload-and-process"
data = {
    "files": ["117_ideal.mol2"],
    "approach": "classical"
}

# Convert to JSON and encode
json_data = json.dumps(data).encode('utf-8')

# Create request
req = urllib.request.Request(url, 
                           data=json_data,
                           headers={'Content-Type': 'application/json'})

try:
    with urllib.request.urlopen(req) as response:
        result = response.read().decode('utf-8')
        print("Response received:")
        print(result[:500])  # First 500 characters
        
        # Parse JSON
        result_data = json.loads(result)
        if 'session_id' in result_data:
            print(f"Session ID: {result_data['session_id']}")
        
except Exception as e:
    print(f"Error: {e}")
