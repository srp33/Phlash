<template>
  <div class="wrapper" @click="checkIfFilesUploaded">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <h1>Upload Files</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p v-if="fasta && genbank">You have uploaded all required files. Click 'Next' to continue.</p>
        <p v-else>Please upload the GenBank file from DNA Master <strong>(.gb, .gbk, .gbf)</strong> and then double click 'Next'</p>
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
          GeneMark is automatically executed once you upload your GenBank file.<br />
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
        <hr />
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
              name: 'Blast',
              params: { phageID: $route.params.phageID },
            }"
            :event="fasta && genbank ? 'click' : ''"
          >
            <button class="btn btn-light btn-nav disabled" id="next-top" @click="checkIfFilesUploaded">
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

      <div class="upload-wrapper genbank">
        <h5 class="upload-title">GenBank file</h5>
        <vue-dropzone ref="myVueDropzone" id="dropzone" :duplicateCheck="true" :options="dropzoneOptions" :destroyDropzone="false"></vue-dropzone>
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
            name: 'Blast',
            params: { phageID: $route.params.phageID },
          }"
          :event="fasta && genbank ? 'click' : ''"
        >
          <button class="btn btn-light btn-nav disabled" id="next-bottom" @click="checkIfFilesUploaded">
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
import vue2Dropzone from "vue2-dropzone";
import 'vue2-dropzone/dist/vue2Dropzone.min.css';

export default {
  name: "Upload",
  components: {
    VueDropzone: vue2Dropzone,
    Loading,
    Navbar,
  },

  data() {

    return {
      fasta: false,
      genbank: false,
      showDelayMessage: true,
      dropzoneOptions: this.setDropzone(),
      blastCompleted: false,
    };

  },

  created() {
    this.checkIfFilesUploaded();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    // navDNAMaster: function () {
    //   if (this.fasta && this.genbank) return true;
    //   else return false;
    // },

    navBlast: function () {
      if (this.fasta && this.genbank) return true;
      else return false;
    },

    navGeneMap: function () {
      if (this.fasta && this.genbank) return true;
      else return false;
    },

    navAnnotations: function () {
      return this.blastCompleted;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
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

    /**
     * Sets all functionality for dropzone.
     * See dropzone specs online for descriptions of the parameters.
     */
    setDropzone() {
      return {
        url: this.getUploadUrl(),
        addRemoveLinks: true,
        acceptedFiles: ".gb, .gbk, .gbf",
        chunking: false,
        maxFiles: 1,
        dictDefaultMessage: "Drag files here or click to browse",
        dictInvalidFileType: "Only '.gb', '.gbk', or '.gbf' file types are allowed.",
        dictRemoveFileConfirmation: "Are you sure you want to remove this file? This will remove all progress that you have made on this Phage.",
        dictMaxFilesExceeded: "You can only upload one file.",
        init: function() {

          axios
            .post(
              this.options.url.slice(0,this.options.url.indexOf("uploadGenbank")) + `display/none`
            )
            .then((response) => {
              console.log(response.data);
              var fileName = response.data.genbank_file;
              var fileSize = response.data.genbank_file_size;
              if (fileName != "Not found") {
                this.addCustomFile(
                  // File options
                  {
                      // flag: processing is complete
                      processing: true,
                      // flag: file is accepted (for limiting maxFiles)
                      accepted: true,
                      // name of file on page
                      name: fileName,
                      // image size
                      size: fileSize,
                      // image type
                      type: '.gb',
                      // flag: status upload
                      status: this.SUCCESS,
                      lastModifiedDate: 'unimportant',
                  },
                  // Custom response for event success
                  {
                      status: "success"
                  }
                );
              }
              console.log(this.files);
            })
            .catch((error) => {
              console.log(error);
            });

          this.addCustomFile = function(file, response){
            // Push file to collection
            this.files.push(file);
            // Emulate event to create interface
            this.emit("addedfile", file);
            // Add status processing to file
            this.emit("processing", file);
            // Add status success to file AND RUN EVENT success from responce
            this.emit("success", file, response , false);
            // Add status complete to file
            this.emit("complete", file);
          }

          this.on("removedfile", function(file) {
            console.log(file);
            if (file.processing) {
            axios
              .post(
                this.options.url.slice(0,this.options.url.indexOf("uploadGenbank")) + `delete/${file.name}`
              )
              .then((response) => {
                console.log(response.data);
              })
              .catch((error) => {
                console.log(error);
              });
            }
          });
        }
      };
    },

    /**
     * Returns the upload URL for dropzone.
     */
    getUploadUrl() {
      return `http://daniel.byu.edu:5000/phlash_api/upload/${this.$route.params.phageID}/uploadGenbank/none`;
    },

    /**
     * Checks to see if fasta and genbank files have been uploaded.
     */
    checkIfFilesUploaded() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/check_upload/${this.$route.params.phageID}`
        )
        .then((response) => {
          console.log(response);
          this.blastCompleted = response.data.blast_completed;
          this.fasta = response.data.fasta;
          this.genbank = response.data.genbank;
          if (this.fasta && this.genbank) {
            this.showDelayMessage = false;
          }
        })
        .catch((error) => {
          console.log(error);
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
