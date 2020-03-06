<template>
   <div class="container">
      <h1>Upload Files</h1>
      <div class="alert alert-primary">
         <p><strong>Instructions</strong></p>
         <p>The following files (with the accepted extensions) are required in order to continue:</p>
         <ul>
            <li v-if="fasta">FASTA file containing the full genome <strong>(.fasta, .fna, .faa, .ffn, .frn)</strong></li>
            <li v-if="genbank">GenBank file from DNA Master <strong>(.gb, .gbk)</strong></li>
            <li v-if="gdata">gdata file from GeneMark <strong>(.gdata)</strong></li>
            <li v-if="ldata">ldata file from GeneMark <strong>(.ldata)</strong></li>
            <li v-if="allUploaded">You have uploaded the required files. Thank you!</li>
         </ul>
         <router-link  to="/dnamaster" v-if="allUploaded">
            <button class="btn btn-primary btn-large" id="next-top"><strong>Next</strong></button>
         </router-link>
      </div>

      <div id="success-alert" class="alert alert-success" role="alert" v-if="showSuccessAlert"></div>
      <div id="danger-alert" class="alert alert-danger" role="alert" v-if="showDangerAlert"></div>

      <div class="upload">
         <form id="upload-form" role="form" enctype="multipart/form-data">
            <div class="upload-btn-wrapper">
               <button class="btn btn-upload">Drag files here or<br/>click to browse!</button>
               <input type="file" id="files" ref="files" name="files" multiple v-on:change="handleFilesUpload" class="form-control">
            </div>
            <div class="selected-files" v-if="showSelectedFiles">
               <strong>SELECTED FILES:</strong>
               <ul>
                  <li v-for="file in files" :key="file">
                     {{ file.name }}
                  </li>
               </ul>
            </div>
         </form>
      </div>

      <button class="btn btn-success btn-upload-submit" v-if="showSelectedFiles" @click="upload">Upload</button>

      <!-- <br>
      <button class="btn btn-primary" @click="updateDatabase">Add DNA Master data to database</button>
      <br><br> -->
      <router-link  to="/dnamaster" v-if="allUploaded">
         <button class="btn btn-primary" id="next-bottom"><strong>Next</strong></button>
      </router-link>
   </div>
</template>

<script>
import axios from 'axios';
import Vue from 'vue';

export default {
   data() {
      return {
         files: [],
         showDangerAlert: false,
         showSuccessAlert: false,
         showSelectedFiles: false,
         fasta: true,
         genbank: true,
         gdata: true,
         ldata: true,
         allUploaded: false,
         showDatabase: false,
      }
   },
   methods: {
      elt(type, ...children) {
         let node = document.createElement(type);
         for (let child of children) {
            if (typeof child != "string") node.appendChild(child);
            else node.appendChild(document.createTextNode(child));
         }
         return node;
      },
      handleFilesUpload() {
         this.files = document.querySelector('#upload-form').files.files
         this.showSelectedFiles = true;
         console.log(this.files)
      },
      upload(e) {
         var data = new FormData(document.querySelector('#upload-form'));
         axios.post('http://localhost:5000/api/upload',
            data,
            {
               headers: {
                  'Content-Type': 'multipart/form-data'
               }
            }
         ).then((response) => {
            console.log(response);
            if (response.data.uploaded.length > 0) {
               this.showSuccessAlert = true;
               let successMessage = "<strong>Files uploaded successfully:</strong> ";
               response.data.uploaded.forEach((filename, index) => {
                  if (index === response.data.uploaded.length - 1) {
                     successMessage += filename;
                  } else {
                     successMessage += filename + ", "
                  }
               });
               Vue.nextTick(() => {
                  document.getElementById("success-alert").innerHTML = successMessage;
               });
            }
            
            if (response.data.not_allowed.length > 0) {
               this.showDangerAlert = true;
               let dangerMessage = "<strong>Unacceptable file types:</strong> ";
               response.data.not_allowed.forEach((filename, index) => {
                  if (index === response.data.not_allowed.length - 1) {
                     dangerMessage += filename;
                  } else {
                     dangerMessage += filename + ", "
                  }
               });
               Vue.nextTick(() => {
                  document.getElementById("danger-alert").innerHTML = dangerMessage;
               });
            }

            if (response.data.required.length > 0) {
               let required = response.data.required;
               if (!required.includes("fasta")) this.fasta = false;
               if (!required.includes("genbank")) this.genbank = false;
               if (!required.includes("gdata")) this.gdata = false;
               if (!required.includes("ldata")) this.ldata = false;
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
  padding-top: 40px;
  padding-bottom: 40px;
  width: 100%;
  border-radius: 8px;
  font-size: 20px;
  font-weight: bold;
}

.upload-btn-wrapper input[type=file] {
  font-size: 100px;
  position: absolute;
  left: 0;
  top: 0;
  opacity: 0;
}

.selected-files {
   display: inline-block;
   text-align: left;
   margin: 10px;
}

.selected-files ul {
   margin: 0;
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
