<template>
  <div class="container">
    <h1>Upload Files</h1>
    <div class="alert alert-primary">
      <p><strong>Instructions</strong></p>
      <p>Please upload the following files:</p>
      <ul>
        <li v-if="!fasta">FASTA file containing the full genome* <strong>(.fasta, .fna)</strong></li>
        <li v-if="!genbank">GenBank file from DNA Master <strong>(.gb, .gbk)</strong></li>
        <li v-if="fasta && genbank">
          You have uploaded all required files. Thank you!
        </li>
      </ul>
      <div class="alert alert-info">
        <svg class="bi bi-info-circle-fill" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
          <path fill-rule="evenodd" d="M8 16A8 8 0 108 0a8 8 0 000 16zm.93-9.412l-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM8 5.5a1 1 0 100-2 1 1 0 000 2z" clip-rule="evenodd"/>
        </svg>
        GeneMark is automatically executed once you upload your FASTA file.<br />
        <svg class="bi bi-info-circle-fill" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
          <path fill-rule="evenodd" d="M8 16A8 8 0 108 0a8 8 0 000 16zm.93-9.412l-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM8 5.5a1 1 0 100-2 1 1 0 000 2z" clip-rule="evenodd"/>
        </svg>
        Files are read and parsed as they are uploaded. Keep this in mind if uploading files take a moment.
      </div>
      <router-link :to="{ name: 'DNAMaster', params: {phageID: $route.params.phageID} }"
        v-if="fasta && genbank">
        <button class="btn btn-light" id="next-top">
          <strong>Next</strong>
        </button>
      </router-link>
    </div>

    <div class="upload-wrapper">
      <h5 class="upload-title">Fasta file</h5>
      <div class="alert alert-success" id="fasta-success-alert" role="alert" v-if="showFastaSuccessAlert"></div>
      <div class="alert alert-danger" id="fasta-danger-alert" role="alert" v-if="showFastaDangerAlert"></div>
      <div class="upload">
        <form id="fasta-upload-form" role="form" enctype="multipart/form-data">
          <div class="upload-btn-wrapper">
            <button class="btn btn-upload">
              <loading :active.sync="fastaLoading" :is-full-page="false" :height="80" :width="80"></loading>
              Drag files here or click to browse <br />
              <div class="selected-file" v-if="showFastaFile">
                <strong>Selected file: {{ this.fastaFile.name }}</strong>
              </div>
            </button>
            <input class="form-control" id="file" type="file" ref="file" name="file" v-on:change="handleFileUpload('fasta')" />
          </div>
        </form>
      </div>
      <button class="btn btn-dark btn-upload-submit" v-if="showFastaFile" @click="uploadFile('fasta')" >
        <strong>Upload</strong>
      </button>
    </div>

    <div class="upload-wrapper genbank">
      <h5 class="upload-title">GenBank file</h5>
      <div class="alert alert-success" id="genbank-success-alert" role="alert" v-if="showGenBankSuccessAlert" ></div>
      <div class="alert alert-danger" id="genbank-danger-alert" role="alert" v-if="showGenBankDangerAlert"></div>
      <div class="upload">
        <form id="genbank-upload-form" role="form" enctype="multipart/form-data">
          <div class="upload-btn-wrapper">
            <button class="btn btn-upload">
              <loading :active.sync="genbankLoading" :is-full-page="false" :height="80" :width="80"></loading>
              Drag files here or click to browse <br />
              <div class="selected-file" v-if="showGenBankFile">
                <strong>Selected file: {{ this.genbankFile.name }}</strong>
              </div>
            </button>
            <input class="form-control" id="file" type="file" ref="file" name="file" v-on:change="handleFileUpload('genbank')" />
          </div>
        </form>
      </div>
      <button class="btn btn-dark btn-upload-submit" v-if="showGenBankFile" @click="uploadFile('genbank')">
        <strong>Upload</strong>
      </button>
    </div>
    
    <router-link :to="{ name: 'DNAMaster', params: {phageID: $route.params.phageID} }" v-if="fasta && genbank">
      <button class="btn btn-light" id="next-bottom">
        <strong>Next</strong>
      </button>
    </router-link>
  </div>
</template>

<script>
import axios from "axios";
import Vue from "vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Upload",
  components: {
    Loading
  },
  data() {
    return {
      fasta: false,
      fastaFile: null,
      fastaLoading: false,
      showFastaFile: false,
      showFastaDangerAlert: false,
      showFastaSuccessAlert: false,
      genbank: false,
      genbankFile: null,
      genbankLoading: false,
      showGenBankFile: false,
      showGenBankDangerAlert: false,
      showGenBankSuccessAlert: false,
    };
  },
  methods: {
    handleFileUpload(fileType) {
      if (fileType === "fasta") {
        this.fastaFile = document.querySelector(`#${fileType}-upload-form`).file.files[0];
        this.showFastaFile = true;
      } else if (fileType === "genbank") {
        this.genbankFile = document.querySelector(`#${fileType}-upload-form`).file.files[0];
        this.showGenBankFile = true;
      }
    },
    uploadFile(fileType, e) {
      var data = new FormData();
      if (fileType === "fasta") {
        this.fastaLoading = true;
        data.append("file", this.fastaFile);
      }
      else if (fileType === "genbank") {
        this.genbankLoading = true;
        data.append("file", this.genbankFile);
      }
      data.append("fileType", fileType);
      axios.post(`http://localhost:5000/api/upload/${this.$route.params.phageID}`,
          data,
          {
            headers: {
              "Content-Type": "multipart/form-data"
            }
          }
        )
        .then(response => {
          if (typeof response.data.uploaded !== "undefined") {
            let fileExt = response.data.uploaded.split(".").pop();
            if (fileExt === "fasta" || fileExt === "fna") {
              this.fastaLoading = false;
              this.showFastaSuccessAlert = true;
              this.showFastaDangerAlert = false;
              this.fasta = true;
            } else if (fileExt === "gb" || fileExt === "gbk") {
              this.genbankLoading = false;
              this.showGenBankSuccessAlert = true;
              this.showGenBankDangerAlert = false;
              this.genbank = true;
            }
            let successMessage = `<strong>${response.data.uploaded}</strong> uploaded successfully!`;
            Vue.nextTick(() => {
              document.getElementById(
                `${fileType}-success-alert`
              ).innerHTML = successMessage;
            });
          }

          if (typeof response.data.not_allowed !== "undefined") {
            let fileExt = response.data.not_allowed.split(".").pop();
            console.log(fileExt);
            if (fileType === "fasta") {
              this.fastaLoading = false;
              this.showFastaDangerAlert = true;
              this.showFastaSuccessAlert = false;
            } else if (fileType === "genbank") {
              this.genbankLoading = false;
              this.showGenBankDangerAlert = true;
              this.showGenBankSuccessAlert = false;
            }
            let dangerMessage = `<strong>${fileExt}</strong> is an unacceptable ${fileType} file extension.`;
            Vue.nextTick(() => {
              document.getElementById(
                `${fileType}-danger-alert`
              ).innerHTML = dangerMessage;
            });
        }
      })
      .catch(error => {
        console.log(error);
      });
    }
  }
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

.bi-info-circle-fill {
  margin-right: 3px;
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

.upload-btn-wrapper input[type="file"] {
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
  margin: auto;
  width: 100%;
}

.genbank {
  margin-bottom: 15px;
}
/* ----- Rest of Page ----- */
#next-top, #next-bottom {
  margin: 10px auto;
}
</style>
