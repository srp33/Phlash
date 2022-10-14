# Phlash - client

## Project setup
- Create a google client ID for the application.
    - Visit [Google Cloud Platform](https://console.cloud.google.com/apis/credentials) and click "Create Credentials" and then "OAuth client ID". 
    - Under "Application type" select "Web Application".
    - Enter a name for the application.
    - Under "Authorized JavaScript origins," click "Add URL" and enter the URL that points to the server you are hosting. 
    - Under "Authorized redirect URIs," click "Add URL" and enter the URL that points to the server you are hosting. 
    - Click "CREATE".
    - Copy your client ID.
- Create a file called `.env.production` in the front-end directory to point to the API (back-end) server and to store your google client id. 
    - Within the file add the two following lines with your url and client ID:  
        - VUE_APP_BASE_URL=https://example.url:5000/phlash_api  
        - VUE_APP_API_KEY=my-clientID.apps.googleusercontent.com  
- Run the `build_docker` script to build the Docker image, then execute the `run_docker` script to start the container for the front-end.
    - The server is configured to run on port 5050 although this may be changed in in the `run_docker` script.
