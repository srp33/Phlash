<template>
  <div class="wrapper" @click="clearNotifications()">
    <Navbar
      :upload="navUpload"
      :dnamaster="navDNAMaster"
      :blast="navBlast"
      :annotations="navAnnotations"
    />
    <div class="container">
      <h1>Upload Files</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>Please upload the following files:</p>
        <ul>
          <li v-if="!fasta">
            FASTA file containing the full genome*
            <strong>(.fasta, .fna)</strong>
          </li>
          <li v-if="!genbank">
            GenBank file from DNA Master <strong>(.gb, .gbk)</strong>
          </li>
          <li v-if="fasta && genbank">
            You have uploaded all required files. Click 'Next' to continue.
          </li>
        </ul>
        <div class="alert alert-warning" v-if="showDelayMessage">
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
          GeneMark is automatically executed once you upload your FASTA file.<br />
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
          Files are read and parsed as they are uploaded. Keep this in mind if
          uploading files take a moment.
        </div>
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Home' }">
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
              name: 'DNAMaster',
              params: { phageID: $route.params.phageID },
            }"
            :event="fasta && genbank ? 'click' : ''"
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

      <div class="upload-wrapper">
        <h5 class="upload-title">Fasta file</h5>
        <div
          class="alert alert-warning"
          id="fasta-warning-alert"
          role="alert"
          v-if="showNoMoreFastaFiles"
        ></div>
        <div
          class="alert alert-success"
          id="fasta-success-alert"
          role="alert"
          v-if="showFastaSuccessAlert"
        ></div>
        <div
          class="alert alert-danger"
          id="fasta-danger-alert"
          role="alert"
          v-if="showFastaDangerAlert"
        ></div>
        <div class="upload">
          <form
            id="fasta-upload-form"
            role="form"
            enctype="multipart/form-data"
          >
            <div class="upload-btn-wrapper">
              <button class="btn btn-upload">
                <loading
                  :active.sync="fastaLoading"
                  :is-full-page="false"
                  :height="80"
                  :width="80"
                ></loading>
                Drag files here or click to browse <br />
                <div class="selected-file" v-if="showFastaFile">
                  <strong>Selected file: {{ this.fastaFile.name }}</strong>
                </div>
              </button>
              <input
                class="form-control"
                id="file"
                type="file"
                ref="file"
                name="file"
                v-on:change="handleFileUpload('fasta')"
              />
            </div>
          </form>
        </div>
        <button
          class="btn btn-dark btn-upload-submit"
          v-if="showFastaFile"
          @click="uploadFile('fasta')"
        >
          <strong>Upload</strong>
        </button>
        <div class="file-list" v-if="fasta">
          <button
            class="btn btn-light btn-download"
            @click="downloadFile(fastaFileName)"
          >
            <strong>Download</strong> {{ fastaFileName }}
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
            @click="showFastaModal = true"
          >
            <strong>Delete</strong> {{ fastaFileName }}
            <svg
              width="2em"
              height="2em"
              viewBox="0 0 16 16"
              class="bi bi-trash"
              fill="this.backgroundColor"
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
        </div>
      </div>

      <div class="upload-wrapper genbank">
        <h5 class="upload-title">GenBank file</h5>
        <div
          class="alert alert-warning"
          id="genbank-warning-alert"
          role="alert"
          v-if="showNoMoreGenbankFiles"
        ></div>
        <div
          class="alert alert-success"
          id="genbank-success-alert"
          role="alert"
          v-if="showGenBankSuccessAlert"
        ></div>
        <div
          class="alert alert-danger"
          id="genbank-danger-alert"
          role="alert"
          v-if="showGenBankDangerAlert"
        ></div>
        <div class="upload">
          <form
            id="genbank-upload-form"
            role="form"
            enctype="multipart/form-data"
          >
            <div class="upload-btn-wrapper">
              <button class="btn btn-upload">
                <loading
                  :active.sync="genbankLoading"
                  :is-full-page="false"
                  :height="80"
                  :width="80"
                ></loading>
                Drag files here or click to browse <br />
                <div class="selected-file" v-if="showGenBankFile">
                  <strong>Selected file: {{ this.genbankFile.name }}</strong>
                </div>
              </button>
              <input
                class="form-control"
                id="file"
                type="file"
                ref="file"
                name="file"
                v-on:change="handleFileUpload('genbank')"
              />
            </div>
          </form>
        </div>
        <button
          class="btn btn-dark btn-upload-submit"
          v-if="showGenBankFile"
          @click="uploadFile('genbank')"
        >
          <strong>Upload</strong>
        </button>
        <div class="file-list" v-if="genbank">
          <button
            class="btn btn-light btn-download"
            @click="downloadFile(genbankFileName)"
          >
            <strong>Download</strong> {{ genbankFileName }}
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
            @click="showGenBankModal = true"
          >
            <strong>Delete</strong> {{ genbankFileName }}
            <svg
              width="2em"
              height="2em"
              viewBox="0 0 16 16"
              class="bi bi-trash"
              fill="this.backgroundColor"
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
        </div>
      </div>
      <div class="nav-btns-wrapper">
        <router-link :to="{ name: 'Home' }">
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
            name: 'DNAMaster',
            params: { phageID: $route.params.phageID },
          }"
          :event="fasta && genbank ? 'click' : ''"
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
    <b-modal
      v-model="showGenBankModal"
      ref="GBModal"
      title="Warning!"
      hide-footer
    >
      <p>
        Deleting this file will remove all progress that you have made on this
        phage.
      </p>
      <hr />
      <b-button
        class="mt-3"
        block
        style="margin-top: 0px"
        @click="deleteFile(genbankFileName, 'gb')"
      >
        <strong>Delete GenBank file</strong>
      </b-button>
    </b-modal>
    <b-modal
      v-model="showFastaModal"
      ref="FastaModal"
      title="Warning!"
      hide-footer
    >
      <p>
        Deleting this file will remove all progress that you have made on this
        phage.
      </p>
      <hr />
      <b-button
        class="mt-3"
        block
        style="margin-top: 0px"
        @click="deleteFile(fastaFileName, 'fasta')"
      >
        <strong>Delete Fasta file</strong>
      </b-button>
    </b-modal>
  </div>
</template>

<script>
import axios from "axios";
import Vue from "vue";
import Navbar from "../components/Navbar.vue";
import Loading from "vue-loading-overlay";
import "vue-loading-overlay/dist/vue-loading.css";

export default {
  name: "Upload",
  components: {
    Loading,
    Navbar,
  },

  data() {

    return {
      fasta: false,
      fastaFile: null,
      fastaFileName: null,
      fastaLoading: false,
      showFastaFile: false,
      showFastaDangerAlert: false,
      showFastaSuccessAlert: false,
      genbank: false,
      genbankFile: null,
      genbankFileName: null,
      genbankLoading: false,
      showGenBankFile: false,
      showGenBankDangerAlert: false,
      showGenBankSuccessAlert: false,
      showDelayMessage: false,
      showNoMoreFastaFiles: false,
      showNoMoreGenbankFiles: false,
      showFastaModal: false,
      showGenBankModal: false,
    };

  },

  created() {
    this.checkIfFilesUploaded();
    this.displayFiles();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    navDNAMaster: function () {
      if (this.fasta && this.genbank) return true;
      else return false;
    },

    navBlast: function () {
      return false;
    },

    navAnnotations: function () {
      return false;
    },

  },

  watch: {

    fasta: function () {
      if (this.fasta && this.genbank) {
        document.getElementById("next-top").classList.remove("disabled");
        document.getElementById("next-bottom").classList.remove("disabled");
      }
    },

    genbank: function () {
      if (this.fasta && this.genbank) {
        document.getElementById("next-top").classList.remove("disabled");
        document.getElementById("next-bottom").classList.remove("disabled");
      }
    },

  },

  methods: {
    clearNotifications() {
      this.showFastaDangerAlert = false;
      this.showFastaSuccessAlert = false;
      this.showGenBankDangerAlert = false;
      this.showGenBankSuccessAlert = false;
    },

    checkIfFilesUploaded() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/check_upload/${this.$route.params.phageID}`
        )
        .then((response) => {
          this.fasta = response.data.fasta;
          this.genbank = response.data.genbank;
        })
        .catch((error) => {
          console.log(error);
        });
    },

    handleFileUpload(fileType) {
      if (fileType === "fasta") {
        this.fastaFile = document.querySelector(
          `#${fileType}-upload-form`
        ).file.files[0];
        this.showFastaFile = true;
      } else if (fileType === "genbank") {
        this.genbankFile = document.querySelector(
          `#${fileType}-upload-form`
        ).file.files[0];
        this.showGenBankFile = true;
      }
    },

    uploadFile(fileType, e) {
      var data = new FormData();
      this.showFastaDangerAlert = false;
      this.showFastaSuccessAlert = false;
      this.showGenBankDangerAlert = false;
      this.showGenBankSuccessAlert = false;
      if (fileType === "fasta") {
        if (this.fastaFileName != null) {
          this.showNoMoreFastaFiles = true;
          let warningMessage =
            "You must remove a file before uploading another.";
          Vue.nextTick(() => {
            document.getElementById(
              "fasta-warning-alert"
            ).innerHTML = warningMessage;
          });
          return;
        }
        this.fastaLoading = true;
        data.append("file", this.fastaFile);
      } else if (fileType === "genbank") {
        if (this.genbankFileName != null) {
          this.showNoMoreGenbankFiles = true;
          let warningMessage =
            "You must remove a file before uploading another.";
          Vue.nextTick(() => {
            document.getElementById(
              "genbank-warning-alert"
            ).innerHTML = warningMessage;
          });
          return;
        }
        this.genbankLoading = true;
        data.append("file", this.genbankFile);
      }
      this.showNoMoreFastaFiles = false;
      this.showNoMoreGenbankFiles = false;
      this.showDelayMessage = true;
      data.append("fileType", fileType);
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/upload/${this.$route.params.phageID}/upload/none`,
          data,
          {
            headers: {
              "Content-Type": "multipart/form-data",
            },
          }
        )
        .then((response) => {
          if (typeof response.data.uploaded !== "undefined") {
            let fileExt = response.data.uploaded.split(".").pop();
            if (fileExt === "fasta" || fileExt === "fna") {
              this.fastaFileName = this.fastaFile.name;
              this.fastaLoading = false;
              this.showFastaSuccessAlert = true;
              this.fasta = true;
            } else if (fileExt === "gb" || fileExt === "gbk") {
              this.genbankFileName = this.genbankFile.name;
              this.genbankLoading = false;
              this.showGenBankSuccessAlert = true;
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
            if (fileType === "fasta") {
              this.fastaLoading = false;
              this.showFastaDangerAlert = true;
            } else if (fileType === "genbank") {
              this.genbankLoading = false;
              this.showGenBankDangerAlert = true;
            }
            let dangerMessage = `<strong>${fileExt}</strong> is an unacceptable ${fileType} file extension.`;
            Vue.nextTick(() => {
              document.getElementById(
                `${fileType}-danger-alert`
              ).innerHTML = dangerMessage;
            });
          }
          this.showDelayMessage = false;
        })
        .catch((error) => {
          console.log(error);
        });
    },

    displayFiles() {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/upload/${this.$route.params.phageID}/display/none`
        )
        .then((response) => {
          console.log(response);
          if (response.data.fasta_file != "Not found")
            this.fastaFileName = response.data.fasta_file;
          console.log(this.fastaFileName);
          if (response.data.genbank_file != "Not found")
            this.genbankFileName = response.data.genbank_file;
        })
        .catch((error) => {
          console.log(error);
        });
    },

    deleteFile(fileName, type) {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/upload/${this.$route.params.phageID}/delete/${fileName}`
        )
        .then((response) => {
          if (type == "gb") {
            this.$refs.GBModal.hide();
            this.genbank = false;
            this.genbankFileName = null;
          } else {
            this.$refs.FastaModal.hide();
            this.fasta = false;
            this.fastaFileName = null;
          }
          console.log(response);
        })
        .catch((error) => {
          console.log(error);
        });
    },

    downloadFile(fileName) {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/upload/${this.$route.params.phageID}/download/${fileName}`
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

.bi-info-circle-fill {
  margin-right: 3px;
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

/* ----- Upload Files ----- */
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

.alert-warning {
  margin-top: 20px;
}
</style>
