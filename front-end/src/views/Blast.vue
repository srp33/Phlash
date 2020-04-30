<template>
<div class="container">
   <h1>BLAST</h1>
   <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>Follow the steps below. Continue when all steps have been completed.</p>
      <router-link :to="{ name: 'Annotations', params: {currentUser: $route.params.currentUser} }" v-show="showBlastSuccessAlert">
         <button class="btn btn-light">
            <strong>Next</strong>
         </button>
      </router-link>
   </div>
   <div class="steps">
      <ol>
         <li class="step">
            Download the FASTA file that will be used
            as input for BLAST.<br/>
            <button class="btn btn-light" @click="downloadFile">
               <strong>Download FASTA file</strong>
            </button>
         </li>
         <li class="step">
            Go to BLASTp's website. <i>At the website, make sure to upload
            your file and set the appropriate parameters, as show in the list
            and screenshot below.</i><br/>
            <ul>
               <li><strong>Upload File:</strong> Upload the FASTA file from step 1.</li>
               <li><strong>Database:</strong> Non-redundant protein sequences (nr)</li>
               <li><strong>Algorithm:</strong> blastp (protein-protein BLAST)</li>
            </ul>
            <button class="btn btn-light" @click="goToNCBI">
               <strong>Go to NCBI's BLASTp</strong>
            </button>
            <img id="step-two" src="/images/blast_step2.png"/>
         </li>
         <li class="step">
            In the top left table on the results page, click on
            <i>"Download All."</i> This will show you file formatting options
            for downloading your results. Choose
            <strong>Single-file JSON</strong>. Continue when you have your
            downloaded file ready for upload.<br/>
            <img id="step-three" src="/images/blast_step3.png"/>
         </li>
         <li class="step">
            Upload your BLAST results.
            <div class="upload-wrapper">
               <div id="blast-success-alert" class="alert alert-success" role="alert" v-if="showBlastSuccessAlert"></div>
               <div id="blast-danger-alert" class="alert alert-danger" role="alert" v-if="showBlastDangerAlert"></div>
               <div class="upload">
                  <form id="blast-upload-form" role="form" enctype="multipart/form-data">
                     <div class="upload-btn-wrapper">
                        <button class="btn btn-upload">Drag files here or click to browse<br/>
                           <div class="selected-file" v-if="showBlastFile">
                              <strong>Selected file: {{ this.blastFile.name }}</strong>
                           </div>
                        </button>
                        <input type="file" id="file" ref="file" name="file" v-on:change="handleFileUpload()" class="form-control">
                     </div>
                  </form>
               </div>
               <button class="btn btn-dark btn-upload-submit" v-if="showBlastFile" @click="uploadFile()"><strong>Upload</strong></button>
            </div>
            <router-link :to="{ name: 'Annotations', params: {currentUser: $route.params.currentUser} }" v-show="showBlastSuccessAlert">
               <button class="btn btn-light" v-show="showBlastSuccessAlert">
                  <strong>Next</strong>
               </button>
            </router-link>
         </li>
      </ol>
   </div>
</div>
</template>

<script>
import axios from 'axios';
import Vue from 'vue';

export default {
   data() {
      return {
         fileDownloaded: false,
         clickedNCBI: false,
         blast: false,
         blastFile: null,
         showBlastFile: false,
         showBlastDangerAlert: false,
         showBlastSuccessAlert: false,
      };
   },
   methods: {
      downloadFile() {
         axios.post(`http://localhost:5000/api/blast/${this.$route.params.currentUser}/download`)
         .then(response => {
            let data = response.data;
            const blob = new Blob([data], { type: 'application/fasta' })
            let link = document.createElement('a')
            link.href = window.URL.createObjectURL(blob)
            link.download = `${this.$route.params.currentUser}_blast.fasta`
            link.click()
            this.fileDownloaded = true;
         });
      },
      goToNCBI() {
         window.open("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins", "_blank")
         this.clickedNCBI = true;
      },
      handleFileUpload() {
         this.blastFile = document.querySelector('#blast-upload-form').file.files[0]
         this.showBlastFile = true;
      },
      uploadFile(e) {
         var data = new FormData();
         data.append('file', this.blastFile);
         // data.append('fileType', fileType)
         axios.post(`http://localhost:5000/api/blast/${this.$route.params.currentUser}/upload`,
            data,
            {
               headers: {
                  'Content-Type': 'multipart/form-data'
               }
            }
         ).then(response => {
            console.log(response);
            if (typeof response.data.uploaded !== "undefined") {
               let fileExt = response.data.uploaded.split('.').pop();
               this.showBlastSuccessAlert = true;
               this.showBlastDangerAlert = false;
               this.blast = true;
               let successMessage = `<strong>${response.data.uploaded}</strong> uploaded successfully!`;
               Vue.nextTick(() => {
                  document.getElementById('blast-success-alert').innerHTML = successMessage;
               });
            } else if (typeof response.data.not_allowed !== "undefined") {
               let fileExt = response.data.not_allowed.split('.').pop();
               console.log(fileExt)
               this.showBlastDangerAlert = true;
               this.showBlastSuccessAlert = false;
               let dangerMessage = `<strong>${fileExt}</strong> is an unacceptable JSON file extension.`;
               Vue.nextTick(() => {
                  document.getElementById('blast-danger-alert').innerHTML = dangerMessage;
               });
            }
         })
         .catch(error => {
            console.log(error)
         });
      },
   },
};
</script>

<style scoped>
/* ----- Title and Alerts ----- */
h1 {
   margin: 40px auto;
}

.alert-primary {
   text-align: left;
   margin: 40px auto;
}

/* ----- Steps ----- */
.steps {
   text-align: left;
}

.step {
   margin-bottom: 20px;
}

.steps img {
   max-width: 100%;
   height: auto;
   margin: 10px;
}

.btn-light {
   margin: 15px auto;
}

/* ----- Upload ----- */
.upload-wrapper {
   margin: 30px auto;
}

.upload-btn-wrapper {
  position: relative;
  overflow: hidden;
  display: inline-block;
  width: 100%;
}

.btn-upload {
  border: 3px dashed gray;
  color: gray;
  background-color: white;
  padding-top: 20px;
  padding-bottom: 20px;
  width: 100%;
  border-radius: 8px;
  font-size: 18px;
  font-weight: bold;
}

.upload-btn-wrapper input[type=file] {
  font-size: 100px;
  position: absolute;
  left: 0;
  top: 0;
  opacity: 0;
}

.selected-file {
   display: inline-block;
   margin: 10px;
}

.selected-file strong {
   color: #474747;
   font-size: 16px;
   text-align: center;
}

.selected-file span {
   margin: 0;
   color: #474747;
   font-size: 16px;
   text-align: left;
}

.btn-upload-submit {
   display: block;
   margin:  auto;
   width: 100%;
}
</style>