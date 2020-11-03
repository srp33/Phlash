<template>
  <div class="wrapper" @click="clearNotifications">
    <Navbar
      :upload="navUpload"
      :dnamaster="navDNAMaster"
      :blast="navBlast"
      :annotations="navAnnotations"
    />
    <div class="container">
      <h1>BLAST</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          Follow the steps below. Click 'Next' when all the steps have been
          completed.
        </p>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{
              name: 'DNAMaster',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-light btn-nav">
              <svg
                class="bi bi-arrow-left"
                width="1em"
                height="1em"
                viewBox="0 0 16 16"
                fill="currentColor"
                xmlns="http://www.w3.org/2000/svg"
              >
                <path
                  fill-rule="evenodd"
                  d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z"
                  clip-rule="evenodd"
                />
                <path
                  fill-rule="evenodd"
                  d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z"
                  clip-rule="evenodd"
                />
              </svg>
              <strong>Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Annotations',
              params: { phageID: $route.params.phageID },
            }"
            :event="blastDownloaded && blastUploaded ? 'click' : ''"
          >
            <button class="btn btn-light btn-nav disabled" id="next-top">
              <strong>Next</strong>
              <svg
                class="bi bi-arrow-right"
                width="1em"
                height="1em"
                viewBox="0 0 16 16"
                fill="currentColor"
                xmlns="http://www.w3.org/2000/svg"
              >
                <path
                  fill-rule="evenodd"
                  d="M10.146 4.646a.5.5 0 01.708 0l3 3a.5.5 0 010 .708l-3 3a.5.5 0 01-.708-.708L12.793 8l-2.647-2.646a.5.5 0 010-.708z"
                  clip-rule="evenodd"
                />
                <path
                  fill-rule="evenodd"
                  d="M2 8a.5.5 0 01.5-.5H13a.5.5 0 010 1H2.5A.5.5 0 012 8z"
                  clip-rule="evenodd"
                />
              </svg>
            </button>
          </router-link>
        </div>
      </div>
      <div class="steps">
        <ol>
          <li class="step">
            <a href="#" @click="downloadInputFiles" class="alert-link"
              >Download the zipped FASTA file(s)</a
            >
            that will be used as input for BLAST.
            <loading
              :active.sync="downloadLoading"
              :is-full-page="false"
              :height="20"
              :width="20"
            ></loading>
            <p class="zipfile-tip">
              <svg
                class="bi bi-info-circle-fill"
                width="1em"
                height="1em"
                viewBox="0 0 16 16"
                fill="currentColor"
                xmlns="http://www.w3.org/2000/svg"
              >
                <path
                  fill-rule="evenodd"
                  d="M8 16A8 8 0 108 0a8 8 0 000 16zm.93-9.412l-2.29.287-.082.38.45.083c.294.07.352.176.288.469l-.738 3.468c-.194.897.105 1.319.808 1.319.545 0 1.178-.252 1.465-.598l.088-.416c-.2.176-.492.246-.686.246-.275 0-.375-.193-.304-.533L8.93 6.588zM8 5.5a1 1 0 100-2 1 1 0 000 2z"
                  clip-rule="evenodd"
                />
              </svg>
              If you can't open your zip file, try using
              <a
                target="_blank"
                rel="noopener noreferrer"
                href="https://www.7-zip.org/"
                >7-Zip</a
              >.
            </p>
            <p v-if="downloadLoading">Downloading...</p>
          </li>
          <li class="step">
            <a href="#" @click="goToNCBI" class="alert-link"
              >Go to BLASTp's website.
            </a>
            <i
              >At the website, make sure to upload your file and set the
              appropriate parameters, as show in the list and screenshot
              below.</i
            ><br />
            <div class="alert alert-warning" v-if="clickedNCBI">
              This will take several minutes. If it seems to be taking a long
              time, it is still probably working correctly. For help
              troubleshooting,
              <a href="#" @click="goToHelp" class="alert-link">visit the FAQ</a
              >.
            </div>
            <ul>
              <li>
                <strong>Upload File:</strong> Upload the FASTA file from step 1.
              </li>
              <li>
                <strong>Database:</strong> Non-redundant protein sequences (nr)
              </li>
              <li>
                <strong>Algorithm:</strong> blastp (protein-protein BLAST)
              </li>
            </ul>
            <img id="step-two" src="/phlash/images/blast_step2.png" />
          </li>
          <li class="step">
            In the top left table on the results page, click on
            <i>"Download All."</i> This will show you file formatting options
            for downloading your results. Choose
            <strong>Single-file JSON</strong>. Continue when you have your
            downloaded file ready for upload.<br />
            <img id="step-three" src="/phlash/images/blast_step3.png" />
          </li>
          <li class="step">
            Upload your <strong>{{ numFiles }} single-file JSON</strong> BLAST
            results.
            <div class="upload-wrapper">
              <div
                class="alert alert-warning"
                id="blast-warning-alert"
                role="alert"
                v-if="showNoMoreFiles"
              ></div>
              <div
                class="alert alert-success"
                id="blast-success-alert"
                role="alert"
                v-if="showBlastSuccessAlert"
              ></div>
              <div
                class="alert alert-danger"
                id="blast-danger-alert"
                role="alert"
                v-if="showBlastDangerAlert"
              ></div>
              <div class="upload">
                <form
                  id="blast-upload-form"
                  role="form"
                  enctype="multipart/form-data"
                >
                  <div class="upload-btn-wrapper">
                    <button class="btn btn-upload btn-step">
                      <loading
                        :active.sync="blastLoading"
                        :is-full-page="false"
                        :height="80"
                        :width="80"
                      ></loading>
                      Drag files here or click to browse <br />
                      <div class="selected-files" v-if="showBlastFiles">
                        <strong>Selected file(s):</strong>
                        <ul>
                          <li v-for="file in blastFiles" v-bind:key="file">
                            {{ file.name }}
                          </li>
                        </ul>
                      </div>
                    </button>
                    <input
                      class="form-control"
                      id="files"
                      type="file"
                      ref="files"
                      name="files"
                      multiple
                      v-on:change="handleFilesUpload()"
                    />
                  </div>
                </form>
              </div>
              <button
                class="btn btn-dark btn-upload-submit"
                v-if="showBlastFiles"
                @click="uploadOutputFiles()"
              >
                <strong>Upload</strong>
              </button>
              <div class="file-list">
                <ul>
                  <li v-for="fileName in fileNames" v-bind:key="fileName">
                    <button
                      class="btn btn-light btn-download"
                      @click="downloadOutputFile(fileName)"
                    >
                      <strong>Download</strong> {{ fileName }}
                      <svg
                        width="2em"
                        height="2em"
                        viewBox="0 0 16 16"
                        class="bi bi-download"
                        fill="currentColor"
                        xmlns="http://www.w3.org/2000/svg"
                      >
                        <path
                          fill-rule="evenodd"
                          d="M.5 9.9a.5.5 0 0 1 .5.5v2.5a1 1 0 0 0 1 1h12a1 1 0 0 0 1-1v-2.5a.5.5 0 0 1 1 0v2.5a2 2 0 0 1-2 2H2a2 2 0 0 1-2-2v-2.5a.5.5 0 0 1 .5-.5z"
                        />
                        <path
                          fill-rule="evenodd"
                          d="M7.646 11.854a.5.5 0 0 0 .708 0l3-3a.5.5 0 0 0-.708-.708L8.5 10.293V1.5a.5.5 0 0 0-1 0v8.793L5.354 8.146a.5.5 0 1 0-.708.708l3 3z"
                        />
                      </svg>
                    </button>
                    <button
                      class="btn btn-light btn-trash"
                      @click="deleteOutputFile(fileName)"
                    >
                      <strong>Delete</strong> {{ fileName }}
                      <svg
                        width="2em"
                        height="2em"
                        viewBox="0 0 16 16"
                        class="bi bi-trash"
                        fill="currentColor"
                        xmlns="http://www.w3.org/2000/svg"
                      >
                        <path
                          d="M5.5 5.5A.5.5 0 0 1 6 6v6a.5.5 0 0 1-1 0V6a.5.5 0 0 1 .5-.5zm2.5 0a.5.5 0 0 1 .5.5v6a.5.5 0 0 1-1 0V6a.5.5 0 0 1 .5-.5zm3 .5a.5.5 0 0 0-1 0v6a.5.5 0 0 0 1 0V6z"
                        />
                        <path
                          fill-rule="evenodd"
                          d="M14.5 3a1 1 0 0 1-1 1H13v9a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V4h-.5a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1H6a1 1 0 0 1 1-1h2a1 1 0 0 1 1 1h3.5a1 1 0 0 1 1 1v1zM4.118 4L4 4.059V13a1 1 0 0 0 1 1h6a1 1 0 0 0 1-1V4.059L11.882 4H4.118zM2.5 3V2h11v1h-11z"
                        />
                      </svg>
                    </button>
                  </li>
                </ul>
              </div>
            </div>
          </li>
        </ol>
      </div>
      <div class="nav-btns-wrapper">
        <router-link
          :to="{
            name: 'DNAMaster',
            params: { phageID: $route.params.phageID },
          }"
        >
          <button class="btn btn-light btn-nav">
            <svg
              class="bi bi-arrow-left"
              width="1em"
              height="1em"
              viewBox="0 0 16 16"
              fill="currentColor"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                fill-rule="evenodd"
                d="M5.854 4.646a.5.5 0 010 .708L3.207 8l2.647 2.646a.5.5 0 01-.708.708l-3-3a.5.5 0 010-.708l3-3a.5.5 0 01.708 0z"
                clip-rule="evenodd"
              />
              <path
                fill-rule="evenodd"
                d="M2.5 8a.5.5 0 01.5-.5h10.5a.5.5 0 010 1H3a.5.5 0 01-.5-.5z"
                clip-rule="evenodd"
              />
            </svg>
            <strong>Back</strong>
          </button>
        </router-link>
        <router-link
          :to="{
            name: 'Annotations',
            params: { phageID: $route.params.phageID },
          }"
          :event="blastDownloaded && blastUploaded ? 'click' : ''"
        >
          <button class="btn btn-light btn-nav disabled" id="next-bottom">
            <strong>Next</strong>
            <svg
              class="bi bi-arrow-right"
              width="1em"
              height="1em"
              viewBox="0 0 16 16"
              fill="currentColor"
              xmlns="http://www.w3.org/2000/svg"
            >
              <path
                fill-rule="evenodd"
                d="M10.146 4.646a.5.5 0 01.708 0l3 3a.5.5 0 010 .708l-3 3a.5.5 0 01-.708-.708L12.793 8l-2.647-2.646a.5.5 0 010-.708z"
                clip-rule="evenodd"
              />
              <path
                fill-rule="evenodd"
                d="M2 8a.5.5 0 01.5-.5H13a.5.5 0 010 1H2.5A.5.5 0 012 8z"
                clip-rule="evenodd"
              />
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
    Navbar,
  },

  data() {

    return {
      downloadLoading: false,
      blastLoading: false,
      clickedNCBI: false,
      blastDownloaded: false,
      blastUploaded: false,
      blastFiles: [],
      blastFile: null,
      showBlastFiles: false,
      showBlastDangerAlert: false,
      showBlastSuccessAlert: false,
      numFiles: null,
      fileNames: [],
    };

  },

  created() {
    this.setNumFiles();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    navDNAMaster: function () {
      return true;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      if (this.blastDownloaded && this.blastUploaded) return true;
      else return false;
    },

  },

  watch: {

    blastDownloaded: function () {
      if (this.blastDownloaded && this.blastUploaded) {
        document.getElementById("next-top").classList.remove("disabled");
        document.getElementById("next-bottom").classList.remove("disabled");
      } else {
        document.getElementById("next-top").classList.add("disabled");
        document.getElementById("next-bottom").classList.add("disabled");
      }
    },

    blastUploaded: function () {
      if (this.blastDownloaded && this.blastUploaded) {
        document.getElementById("next-top").classList.remove("disabled");
        document.getElementById("next-bottom").classList.remove("disabled");
      } else {
        document.getElementById("next-top").classList.add("disabled");
        document.getElementById("next-bottom").classList.add("disabled");
      }
    },

  },

  methods: {

    /**
     * Sets all of the notification variables to false.
     */
    clearNotifications() {
      this.showBlastDangerAlert = false;
      this.showBlastSuccessAlert = false;
      this.showNoMoreFiles = false;
    },

    /**
     * Checks to see if the blast output files have been uploaded.
     * Checks to see if the blast input folder has been downloaded.
     */
    checkFiles() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/none/none`
        )
        .then((response) => {
          this.blastDownloaded = response.data.blast_downloaded;
          if (this.fileNames.length == this.numFiles) this.blastUploaded = true;
          if (!this.blastDownloaded) this.downloadInputFiles();
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Downloads the blast input files.
     */
    downloadInputFiles() {
      this.downloadLoading = true;
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/downloadInput/none`,
          FormData,
          {
            responseType: "blob",
          }
        )
        .then((response) => {
          console.log(response.data);
          let file_data = response.data;
          const blob = new Blob([file_data], { type: "application/zip" });
          let link = document.createElement("a");
          link.href = window.URL.createObjectURL(blob);
          link.download = `${this.$route.params.phageID}_blast.zip`;
          link.click();
          this.downloadLoading = false;
          this.blastDownloaded = true;
          this.setNumFiles();
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Sets the number of blast output files that must be uploaded.
     */
    setNumFiles() {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/numFiles/none`
        )
        .then((response) => {
          this.numFiles = Number(response.data);
          this.displayOutputFiles();
        });
    },

    /**
     * Links to the NCBI BLASTp page.
     */
    goToNCBI() {
      window.open(
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins",
        "_blank"
      );
      this.clickedNCBI = true;
    },

    /**
     * Links to the NCBI BLAST FAQ page.
     */
    goToHelp() {
      window.open(
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs",
        "_blank"
      );
    },

    /**
     * Gets all of the names of all of the blast output files that have been uploaded.
     */
    displayOutputFiles() {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/displayOutput/none`
        )
        .then((response) => {
          console.log(response);
          this.fileNames = response.data.file_names;
          this.checkFiles();
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Deletes a given file.
     * @param {string} fileName the name of the file to be deleted.
     */
    deleteOutputFile(fileName) {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/deleteOutput/${fileName}`
        )
        .then((response) => {
          var index = this.fileNames.indexOf(fileName);
          this.fileNames.splice(index, 1);
          if (this.fileNames != this.numFiles) this.blastUploaded = false;
          console.log(response);
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Adds the file so that it can be uploaded.
     */
    handleFilesUpload() {
      this.blastFile = document.querySelector("#blast-upload-form").files.files;
      for (var i = 0; i < this.blastFile.length; i++) {
        this.blastFiles.push(this.blastFile[i]);
      }
      this.showBlastFiles = true;
      console.log(this.blastFiles);
    },

    /**
     * Downloads a given file.
     * @param {string} fileName the name of the file to be downloaded.
     */
    downloadOutputFile(fileName) {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/downloadOutput/${fileName}`
        )
        .then((response) => {
          let data = response.data.file_data;
          const blob = new Blob([data], { type: "application/fasta" });
          let link = document.createElement("a");
          link.href = window.URL.createObjectURL(blob);
          link.download = `${fileName}`;
          link.click();
          this.downloadLoading = false;
          this.fileDownloaded = true;
        });
    },

    /**
     * Uploads the blast output files.
     * Handles incorrect file uploads.
     */
    uploadOutputFiles() {
      this.blastLoading = true;
      var data = new FormData();
      var uploadFile = this.blastFiles[this.blastFiles.length - 1];
      if (this.fileNames.length >= this.numFiles) {
        this.blastLoading = false;
        this.showBlastDangerAlert = false;
        this.showBlastSuccessAlert = false;
        this.showNoMoreFiles = true;
        let warningMessage = "You must remove a file before uploading another.";
        Vue.nextTick(() => {
          document.getElementById(
            "blast-warning-alert"
          ).innerHTML = warningMessage;
        });
        return;
      }
      for (var i = 0; i < this.fileNames.length; i++) {
        if (uploadFile.name == this.fileNames[i]) {
          this.blastLoading = false;
          this.showBlastDangerAlert = false;
          this.showBlastSuccessAlert = false;
          this.showNoMoreFiles = true;
          this.blastFiles.splice(this.blastFiles.length - 1, 1);
          let warningMessage = "This file has already been uploaded.";
          Vue.nextTick(() => {
            document.getElementById(
              "blast-warning-alert"
            ).innerHTML = warningMessage;
          });
          return;
        }
      }
      data.append("file", uploadFile);
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/uploadOutput/none`,
          data,
          {
            headers: {
              "Content-Type": "multipart/form-data",
            },
          }
        )
        .then((response) => {
          console.log(response);
          if (typeof response.data.uploaded !== "undefined") {
            this.fileNames.push(uploadFile.name);
            this.blastLoading = false;
            this.showBlastSuccessAlert = true;
            this.showBlastDangerAlert = false;
            this.blast = true;
            if (this.fileNames.length == this.numFiles)
              this.blastUploaded = true;
            let successMessage = `<strong>${response.data.uploaded}</strong> uploaded successfully!`;
            Vue.nextTick(() => {
              document.getElementById(
                "blast-success-alert"
              ).innerHTML = successMessage;
            });
          } else if (typeof response.data.not_allowed !== "undefined") {
            let fileExt = response.data.not_allowed.split(".").pop();
            this.blastLoading = false;
            this.showBlastDangerAlert = true;
            this.showBlastSuccessAlert = false;
            let dangerMessage = `<strong>${fileExt}</strong> is an unacceptable JSON file extension.`;
            Vue.nextTick(() => {
              document.getElementById(
                "blast-danger-alert"
              ).innerHTML = dangerMessage;
            });
          }
        })
        .catch((error) => {
          console.log(error);
        });
      this.blastFiles.splice(this.blastFiles.length - 1, 1);
      if (this.blastFiles.length > 0) {
        this.uploadOutputFiles();
      }
    },

  },
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

.alert-warning {
  margin-top: 20px;
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
.selected-files ul {
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

.btn-download:hover {
  color: rgb(26, 87, 26);
}

.btn-trash:hover {
  color: rgb(143, 27, 27);
}

.btn-download {
  width: 48%;
  margin: 5px;
}

.btn-trash {
  width: 48%;
  margin: 5px;
}

.file-list {
  padding-top: 20px;
}
</style>
