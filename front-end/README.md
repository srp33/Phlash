# Phlash - client

## Project setup
- Create a google client ID for the application.
    - Visit [Google Cloud Platform](https://console.cloud.google.com/apis/credentials) and click "Create Credentials" and then "OAuth client ID". 
    - Enter the URL that points to the server under "Authorized JavaScript origins". 
    - Copy your client ID.
- Create a file called `.env.production` in the front-end directory to point to the API (back-end) server you are hosting and to store your google client id. 
    - Within the file add the two following lines with your url and client ID:  
        - VUE_APP_BASE_URL=https://example.url:5000/phlash_api  
        - GOOGLE=my-clientID.apps.googleusercontent.com  
- Run the `build_docker` script to build the Docker image, then execute the `run_docker` script to start the container for the front-end.
    - The server is configured to run on port 5050 although this may be changed in in the `run_docker` script.
