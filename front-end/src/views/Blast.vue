<template>
  <div class="wrapper">
    <Navbar :upload="navUpload" :dnamaster="navDNAMaster" :blast="navBlast" :annotations="navAnnotations" />
    <div class="container">
      <h1>BLAST</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>Follow the steps below. Continue when all steps have been completed.</p>
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'DNAMaster', params: {phageID: $route.params.phageID} }">
            <button class="btn btn-light btn-nav">
              <svg class="bi bi-arrow-left" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z" clip-rule="evenodd"/>
                <path fill-rule="evenodd" d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z" clip-rule="evenodd"/>
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
          <router-link :to="{ name: 'Annotations', params: {phageID: $route.params.phageID} }"
                       :event="(blastDownloaded && blastUploaded) ? 'click' : ''">
            <button class="btn btn-light btn-nav disabled" id="next-top">
              <strong>Next</strong>
              <svg class="bi bi-arrow-right" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M10.146 4.646a.5.5 0 01.708 0l3 3a.5.5 0 010 .708l-3 3a.5.5 0 01-.708-.708L12.793 8l-2.647-2.646a.5.5 0 010-.708z" clip-rule="evenodd"/>
                <path fill-rule="evenodd" d="M2 8a.5.5 0 01.5-.5H13a.5.5 0 010 1H2.5A.5.5 0 012 8z" clip-rule="evenodd"/>
              </svg>
            </button>
          </router-link>
        </div>
      </div>
      <div class="steps">
        <ol>
          <li class="step">
            Download the zipped FASTA file(s) that will be used as input for BLAST.
            <p class="zipfile-tip">
              <svg class="bi bi-info-circle-fill" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
                <path fill-rule="evenodd" d="M8 16A8 8 0 108 0a8 8 0 000 16zm.93-9.412l-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM8 5.5a1 1 0 100-2 1 1 0 000 2z" clip-rule="evenodd"/>
              </svg> If you can't open your zip file, try using <a target="_blank" rel="noopener noreferrer" href="https://www.7-zip.org/">7-Zip</a>.
            </p>
            <button class="btn btn-light btn-step" style="position: relative;" @click="downloadFile">
              <loading :active.sync="downloadLoading" :is-full-page="false" :height="20" :width="20"></loading>
              <strong>Download FASTA file</strong>
            </button>
            <p v-if="downloadLoading">Downloading...</p>
          </li>
          <li class="step">
            Go to BLASTp's website. <i>At the website, make sure to upload
              your file and set the appropriate parameters, as show in the list
              and screenshot below.</i><br />
            <ul>
              <li><strong>Upload File:</strong> Upload the FASTA file from step 1.</li>
              <li><strong>Database:</strong> Non-redundant protein sequences (nr)</li>
              <li><strong>Algorithm:</strong> blastp (protein-protein BLAST)</li>
            </ul>
            <button class="btn btn-light btn-step" @click="goToNCBI">
              <strong>Go to NCBI's BLASTp</strong>
            </button>
            <!--<img id="step-two" src="<%= BASE_URL %>/images/blast_step2.png" />-->
            <img id="step-two" src="/phlash/images/blast_step2.png" />
          </li>
          <li class="step">
            In the top left table on the results page, click on
            <i>"Download All."</i> This will show you file formatting options for
            downloading your results. Choose <strong>Single-file JSON</strong>. 
            Continue when you have your downloaded file ready for upload.<br />
            <!--<img id="step-three" src="<%= BASE_URL %>/images/blast_step3.png" />-->
            <img id="step-three" src="/phlash/images/blast_step3.png" />
          </li>
          <li class="step">
            Upload your <strong>{{ numBlastFilesDownloaded }} single-file JSON</strong> BLAST results.
            <div class="upload-wrapper">
              <div class="alert alert-success" id="blast-success-alert" role="alert" v-if="showBlastSuccessAlert"></div>
              <div class="alert alert-danger" id="blast-danger-alert" role="alert" v-if="showBlastDangerAlert"></div>
              <div class="upload">
                <form id="blast-upload-form" role="form" enctype="multipart/form-data">
                  <div class="upload-btn-wrapper">
                    <button class="btn btn-upload btn-step">
                      <loading :active.sync="blastLoading" :is-full-page="false" :height="80" :width="80"></loading>
                      Drag files here or click to browse <br />
                      <div class="selected-files" v-if="showBlastFiles">
                        <strong>Selected file(s):</strong>
                        <ul>
                          <li v-for="file in blastFiles" :key="file.name">
                            {{ file.name }}
                          </li>
                        </ul>
                      </div>
                    </button>
                    <input class="form-control" id="files" type="file" ref="files" name="files" multiple v-on:change="handleFilesUpload()" />
                  </div>
                </form>
              </div>
              <button class="btn btn-dark btn-upload-submit" v-if="showBlastFiles" @click="uploadFiles()">
                <strong>Upload</strong>
              </button>
            </div>
          </li>
        </ol>
      </div>
      <div class="nav-btns-wrapper">
        <router-link :to="{ name: 'DNAMaster', params: {phageID: $route.params.phageID} }">
          <button class="btn btn-light btn-nav">
            <svg class="bi bi-arrow-left" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
              <path fill-rule="evenodd" d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z" clip-rule="evenodd"/>
              <path fill-rule="evenodd" d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z" clip-rule="evenodd"/>
            </svg>
            <strong>Back</strong>
          </button>
        </router-link>
        <router-link :to="{ name: 'Annotations', params: {phageID: $route.params.phageID} }"
                      :event="(blastDownloaded && blastUploaded) ? 'click' : ''">
          <button class="btn btn-light btn-nav disabled" id="next-bottom">
            <strong>Next</strong>
            <svg class="bi bi-arrow-right" width="1em" height="1em" viewBox="0 0 16 16" fill="currentColor" xmlns="http://www.w3.org/2000/svg">
              <path fill-rule="evenodd" d="M10.146 4.646a.5.5 0 01.708 0l3 3a.5.5 0 010 .708l-3 3a.5.5 0 01-.708-.708L12.793 8l-2.647-2.646a.5.5 0 010-.708z" clip-rule="evenodd"/>
              <path fill-rule="evenodd" d="M2 8a.5.5 0 01.5-.5H13a.5.5 0 010 1H2.5A.5.5 0 012 8z" clip-rule="evenodd"/>
            </svg>
          </button>
        </router-link>
      </div>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Vue from "vue";
import Navbar from "../components/Navbar.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Blast",
  components: {
    Loading,
    Navbar
  },
  data() {
    return {
      downloadLoading: false,
      blastLoading: false,
      clickedNCBI: false,
      blastDownloaded: false,
      blastUploaded: false,
      blastFiles: null,
      numBlastFilesDownloaded: null,
      showBlastFiles: false,
      showBlastDangerAlert: false,
      showBlastSuccessAlert: false
    };
  },
  created() {
    this.checkFiles();
  },
  computed: {
    navUpload: function() {
      return true;
    },
    navDNAMaster: function() {
      return true;
    },
    navBlast: function() {
      return true;
    },
    navAnnotations: function() {
      if (this.blastDownloaded && this.blastUploaded) return true;
      else return false;
    },
  },
  watch: {
    blastDownloaded: function() {
      if (this.blastDownloaded && this.blastUploaded) {
        document.getElementById("next-top").classList.remove("disabled");
        document.getElementById("next-bottom").classList.remove("disabled");
      }
    },
    blastUploaded: function() {
      if (this.blastDownloaded && this.blastUploaded) {
        document.getElementById("next-top").classList.remove("disabled");
        document.getElementById("next-bottom").classList.remove("disabled");
      }
    },
  },
  methods: {
    checkIfFileUploaded() {
      axios.get(process.env.VUE_APP_BASE_URL + `/blast/${this.$route.params.phageID}/none`)
        .then(response => {
          this.blastDownloaded = response.data.blast_downloaded;
          this.blastUploaded = response.data.blast_uploaded;
          if (!this.blastDownloaded) this.downloadFile();
        })
        .catch(error => {
          console.log(error);
        })
    },
    downloadFile() {
      this.downloadLoading = true;
      axios.post(process.env.VUE_APP_BASE_URL + `/blast/${this.$route.params.phageID}/download`)
        .then(response => {
          console.log(response.data)
          // this.numBlastFilesDownloaded = response.data[0];
          let file_data = response.data;
          const blob = new Blob([file_data], { type: "application/zip" });
          let link = document.createElement("a");
          link.href = window.URL.createObjectURL(blob);
          link.download = `${this.$route.params.phageID}_blast.zip`;
          link.click();
          this.downloadLoading = false;
          this.blastDownloaded = true;
        })
        .catch(error => {
          console.log(error)
        });
    },
    goToNCBI() {
      window.open(
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins",
        "_blank"
      );
      this.clickedNCBI = true;
    },
    handleFilesUpload() {
      this.blastFiles = document.querySelector('#blast-upload-form').files.files
      this.showBlastFiles = true;
    },
    uploadFiles(e) {
      this.blastLoading = true;
      var data = new FormData();
      data.append("file", this.blastFile);
      axios.post(process.env.VUE_APP_BASE_URL + `/blast/${this.$route.params.phageID}/upload`,
          data,
          {
            headers: {
              "Content-Type": "multipart/form-data"
            }
          }
        )
        .then(response => {
          console.log(response);
          if (response.data.uploaded.length > 0) {
            this.blastLoading = false;
            this.showBlastSuccessAlert = true;
            this.blastUploaded = true;
            let successMessage = `<strong>Files uploaded successfully: </strong>`;
            response.data.uploaded.forEach((filename, index) => {
              if (index === response.data.uploaded.length - 1) successMessage += filename;
              else successMessage += filename + ", ";
            });
            Vue.nextTick(() => {
              document.getElementById("blast-success-alert").innerHTML = successMessage;
            });
          }
          
          if (response.data.not_allowed.length > 0) {
            this.blastLoading = false;
            this.showBlastDangerAlert = true;
            let dangerMessage = `<strong>Unacceptable files: </strong>`;
            response.data.not_allowed.forEach((filename, index) => {
              if (index === response.data.not_allowed.length - 1) dangerMessage += filename;
              else dangerMessage += filename + ", ";
            });
            Vue.nextTick(() => {
              document.getElementById("blast-danger-alert").innerHTML = dangerMessage;
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
.wrapper {
  margin: 0;
}

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

.zipfile-tip {
  font-size: 0.9em;
  margin-bottom: 0;
}

.btn-step {
  margin: 15px auto;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 5px;
}

.bi-arrow-left {
  margin-right: 5px;
  margin-left: 0px;
}

.bi-arrow-right {
  margin-right: 0px;
  margin-left: 5px;
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

.upload-btn-wrapper input[type="file"] {
  font-size: 100px;
  position: absolute;
  left: 0;
  top: 0;
  opacity: 0;
}

.selected-files {
  display: inline-block;
  margin: 10px;
}

.selected-files strong, 
.selected-files ul{
  color: #474747;
  font-size: 16px;
  text-align: center;
}

.selected-files ul {
  padding-left: 0;
}

.selected-files span {
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
</style>
