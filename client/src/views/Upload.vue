<template>
<div class="container">
   <h1>Upload Files</h1>
   <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>Please upload the following files:</p>
      <ul>
         <li v-if="!fasta">FASTA file containing the full genome <strong>(.fasta, .fna)</strong></li>
         <li v-if="!genbank">GenBank file from DNA Master <strong>(.gb, .gbk)</strong></li>
         <li v-if="!gdata">gdata file from GeneMark <strong>(.gdata)</strong></li>
         <li v-if="!ldata">ldata file from GeneMark <strong>(.ldata)</strong></li>
         <li v-if="fasta && genbank && gdata && ldata">You have uploaded all required files. Thank you!</li>
      </ul>
      <p>Uploaded files:</p>
      <ul>
         <li v-if="fasta">FASTA file containing the full genome <strong>(.fasta, .fna)</strong></li>
         <li v-if="genbank">GenBank file from DNA Master <strong>(.gb, .gbk)</strong></li>
         <li v-if="gdata">gdata file from GeneMark <strong>(.gdata)</strong></li>
         <li v-if="ldata">ldata file from GeneMark <strong>(.ldata)</strong></li>
         <li v-if="!fasta && !genbank && !gdata && !ldata">You have not uploaded any files yet.</li>
      </ul>
      <router-link  :to="{ name: 'DNAMaster', params: {currentUser: $route.params.currentUser} }" v-if="fasta && genbank && gdata && ldata">
         <button class="btn btn-light" id="next-top"><strong>Next</strong></button>
      </router-link>
   </div>

   <div class="upload-wrapper">
      <h5 class="upload-title">Fasta file</h5>
      <div id="fasta-success-alert" class="alert alert-success" role="alert" v-if="showFastaSuccessAlert"></div>
      <div id="fasta-danger-alert" class="alert alert-danger" role="alert" v-if="showFastaDangerAlert"></div>
      <div class="upload">
         <form id="fasta-upload-form" role="form" enctype="multipart/form-data">
            <div class="upload-btn-wrapper">
               <button class="btn btn-upload">Drag files here or click to browse<br/>
                  <div class="selected-file" v-if="showFastaFile">
                     <strong>Selected file: {{ this.fastaFile.name }}</strong>
                  </div>
               </button>
               <input type="file" id="file" ref="file" name="file" v-on:change="handleFileUpload('fasta')" class="form-control">
            </div>
         </form>
      </div>
      <button class="btn btn-dark btn-upload-submit" v-if="showFastaFile" @click="uploadFile('fasta')"><strong>Upload</strong></button>
   </div>

   <div class="upload-wrapper">
      <h5 class="upload-title">GenBank file</h5>
      <div id="genbank-success-alert" class="alert alert-success" role="alert" v-if="showGenBankSuccessAlert"></div>
      <div id="genbank-danger-alert" class="alert alert-danger" role="alert" v-if="showGenBankDangerAlert"></div>
      <div class="upload">
         <form id="genbank-upload-form" role="form" enctype="multipart/form-data">
            <div class="upload-btn-wrapper">
               <button class="btn btn-upload">Drag files here or click to browse<br/>
                  <div class="selected-file" v-if="showGenBankFile">
                     <strong>Selected file: {{ this.genbankFile.name }}</strong>
                  </div>
               </button>
               <input type="file" id="file" ref="file" name="file" v-on:change="handleFileUpload('genbank')" class="form-control">
            </div>
         </form>
      </div>
      <button class="btn btn-dark btn-upload-submit" v-if="showGenBankFile" @click="uploadFile('genbank')"><strong>Upload</strong></button>
   </div>

   <div class="upload-wrapper">
      <h5 class="upload-title">GData file</h5>
      <div id="gdata-success-alert" class="alert alert-success" role="alert" v-if="showGdataSuccessAlert"></div>
      <div id="gdata-danger-alert" class="alert alert-danger" role="alert" v-if="showGdataDangerAlert"></div>
      <div class="upload">
         <form id="gdata-upload-form" role="form" enctype="multipart/form-data">
            <div class="upload-btn-wrapper">
               <button class="btn btn-upload">Drag files here or click to browse<br/>
                  <div class="selected-file" v-if="showGdataFile">
                     <strong>Selected file: {{ this.gdataFile.name }}</strong>
                  </div>
               </button>
               <input type="file" id="file" ref="file" name="file" v-on:change="handleFileUpload('gdata')" class="form-control">
            </div>
         </form>
      </div>
      <button class="btn btn-dark btn-upload-submit" v-if="showGdataFile" @click="uploadFile('gdata')"><strong>Upload</strong></button>
   </div>

   <div class="upload-wrapper">
      <h5 class="upload-title">LData file</h5>
      <div id="ldata-success-alert" class="alert alert-success" role="alert" v-if="showLdataSuccessAlert"></div>
      <div id="ldata-danger-alert" class="alert alert-danger" role="alert" v-if="showLdataDangerAlert"></div>
      <div class="upload">
         <form id="ldata-upload-form" role="form" enctype="multipart/form-data">
            <div class="upload-btn-wrapper">
               <button class="btn btn-upload">Drag files here or click to browse<br/>
                  <div class="selected-file" v-if="showLdataFile">
                     <strong>Selected file: {{ this.ldataFile.name }}</strong>
                  </div>
               </button>
               <input type="file" id="file" ref="file" name="file" v-on:change="handleFileUpload('ldata')" class="form-control">
            </div>
         </form>
      </div>
      <button class="btn btn-dark btn-upload-submit" v-if="showLdataFile" @click="uploadFile('ldata')"><strong>Upload</strong></button>
   </div>
   <router-link  :to="{ name: 'DNAMaster', params: {currentUser: $route.params.currentUser} }" v-if="fasta && genbank && gdata && ldata">
      <button class="btn btn-light" id="next-top"><strong>Next</strong></button>
   </router-link>
</div>
</template>

<script>
import axios from 'axios';
import Vue from 'vue';

export default {
   data() {
      return {
         fasta: false,
         fastaFile: null,
         showFastaFile: false,
         showFastaDangerAlert: false,
         showFastaSuccessAlert: false,
         genbank: false,
         genbankFile: null,
         showGenBankFile: false,
         showGenBankDangerAlert: false,
         showGenBankSuccessAlert: false,
         gdata: false,
         gdataFile: null,
         showGdataFile: false,
         showGdataDangerAlert: false,
         showGdataSuccessAlert: false,
         ldata: false,
         ldataFile: null,
         showLdataFile: false,
         showLdataDangerAlert: false,
         showLdataSuccessAlert: false,
      }
   },
   methods: {
      handleFileUpload(fileType) {
         if (fileType === "fasta") {
            this.fastaFile = document.querySelector(`#${fileType}-upload-form`).file.files[0]
            this.showFastaFile = true;
         } else if (fileType === "genbank") {
            this.genbankFile = document.querySelector(`#${fileType}-upload-form`).file.files[0]
            this.showGenBankFile = true;
         } else if (fileType === "gdata") {
            this.gdataFile = document.querySelector(`#${fileType}-upload-form`).file.files[0]
            this.showGdataFile = true;
         } else if (fileType === "ldata") {
            this.ldataFile = document.querySelector(`#${fileType}-upload-form`).file.files[0]
            this.showLdataFile = true;
         }
      },
      uploadFile(fileType, e) {
         var data = new FormData();
         if (fileType === "fasta") data.append('file', this.fastaFile);
         else if (fileType === "genbank") data.append('file', this.genbankFile);
         else if (fileType === "gdata") data.append('file', this.gdataFile);
         else if (fileType === "ldata") data.append('file', this.ldataFile);
         data.append('fileType', fileType)
         axios.post(`http://localhost:5000/api/upload/${this.$route.params.currentUser}`,
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
               if (fileExt === "fasta" || fileExt === "fna") {
                  this.showFastaSuccessAlert = true;
                  this.showFastaDangerAlert = false;
                  this.fasta = true;
               } else if (fileExt === "gb" || fileExt === "gbk") {
                  this.showGenBankSuccessAlert = true;
                  this.showGenBankDangerAlert = false;
                  this.genbank = true;
               } else if (fileExt === "gdata") {
                  this.showGdataSuccessAlert = true;
                  this.showGdataDangerAlert = false;
                  this.gdata = true;
               } else if (fileExt === "ldata") {
                  this.showLdataSuccessAlert = true;
                  this.showLdataDangerAlert = false;
                  this.ldata = true;
               }
               let successMessage = `<strong>${response.data.uploaded}</strong> uploaded successfully!`;
               Vue.nextTick(() => {
                  document.getElementById(`${fileType}-success-alert`).innerHTML = successMessage;
               });
            }
            
            if (typeof response.data.not_allowed !== "undefined") {
               let fileExt = response.data.not_allowed.split('.').pop();
               console.log(fileExt)
               if (fileType === "fasta") {
                  this.showFastaDangerAlert = true;
                  this.showFastaSuccessAlert = false;
               } else if (fileType === "genbank") {
                  this.showGenBankDangerAlert = true;
                  this.showGenBankSuccessAlert = false;
               } else if (fileType === "gdata") {
                  this.showGdataDangerAlert = true;
                  this.showGdataSuccessAlert = false;
               } else if (fileType === "ldata") {
                  this.showLdataDangerAlert = true;
                  this.showLdataSuccessAlert = false;
               }
               let dangerMessage = `<strong>${fileExt}</strong> is an unacceptable ${fileType} file extension.`;
               Vue.nextTick(() => {
                  document.getElementById(`${fileType}-danger-alert`).innerHTML = dangerMessage;
               });
            }
            if (response.data.required.length > 0) {
               let required = response.data.required;
               if (!required.includes("fasta")) this.fasta = false;
               if (!required.includes("genbank")) this.genbank = false;
               if (!required.includes("gdata")) this.gdata = false;
               if (!required.includes("ldata")) this.ldata = false;
               if (!this.fasta && !this.genbankFile && !this.gdata && !this.ldata) {
                  this.allUploaded = true;
               }
            } else {
               this.allUploaded = true;
               this.fasta = false;
               this.genbank = false;
               this.gdata = false;
               this.ldata = false;
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

/* ----- Upload ----- */
.upload-wrapper {
   margin: 50px auto;
}

.upload-title {
   margin: 15px auto;
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
   /* text-align: left; */
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

/* ----- Rest of Page ----- */
#next-top {
   margin: 10px auto;
}

#next-bottom {
   margin: 40px auto;
}
</style>
