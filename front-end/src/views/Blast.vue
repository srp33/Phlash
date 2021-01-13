<template>
  <div class="wrapper" @click="checkFiles">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <h1>BLAST</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p>
          Follow the steps below. Double click 'Next' when all the steps have been
          completed.
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{
              name: 'Upload',
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
            <button class="btn btn-light btn-nav disabled" id="next-top"  @click="displayOutputFiles()">
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
        <loading
          :active.sync="downloadLoading"
          :is-full-page="true"
          :height="100"
          :width="100"
        ></loading>
        <ol>
          <li class="step">
            <a href="#" @click="downloadInputFiles" class="alert-link"
              >Download the zipped FASTA file(s)</a
            >
            that will be used as input for BLAST.
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
              >At the website, make sure to upload the files and set the
              appropriate parameters, as show in the list and screenshot
              below. Note that if pop-ups are allowed for this page when 
              the link is clicked the number of tabs that are needed will 
              be opened.</i
            ><br />
            <div class="alert alert-warning alert-dismissible" v-if="clickedNCBI">
              <a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>
              This will take several minutes. If it seems to be taking a long
              time, it is still probably working correctly. <br />
              If an error occurs refresh the page or try re-running that file. Oftentimes it will then work correctly.<br />
              For help troubleshooting,
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
            <vue-dropzone v-if="!blastUploaded" ref="myVueDropzone" id="dropzone" :duplicateCheck="true" :options="dropzoneOptions" :destroyDropzone="false"></vue-dropzone>
            <button v-else class = "btn btn-outline-secondary btn-lg btn-block" style="margin: 30px auto;" @click="removeAll()">BLAST files have been uploaded, click to remove and reupload.</button>
          </li>
        </ol>
      </div>
      <div @mouseenter="displayOutputFiles()" class="nav-btns-wrapper">
        <router-link
          :to="{
            name: 'Upload',
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
          <button class="btn btn-light btn-nav disabled" id="next-bottom" @click="uploadReminder()">
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
  name: "Blast",
  components: {
    VueDropzone: vue2Dropzone,
    Loading,
    Navbar,
  },

  data() {

    return {
      downloadLoading: false,
      clickedNCBI: false,
      blastDownloaded: false,
      blastUploaded: false,
      numFiles: 1,
      fileNames: [],
      dropzoneOptions: this.setDropzone(),
    };

  },

  created() {
    this.setNumFiles();
  },

  computed: {

    navUpload: function () {
      return true;
    },

    // navDNAMaster: function () {
    //   return true;
    // },

    navGeneMap: function () {
      return true;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      if (this.blastDownloaded && this.blastUploaded) return true;
      else return false;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
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
     * Sets all functionality for dropzone.
     * See dropzone specs online for descriptions of the parameters.
     */
    setDropzone() {
      return {
        url: this.getUploadUrl(),
        addRemoveLinks: true,
        acceptedFiles: 'application/json',
        chunking: true,
        maxFiles: null,
        parallelUploads: 5,
        chunkSize: 10000000,
        dictDefaultMessage: "Drag files here or click to browse",
        dictInvalidFileType: "Only '.json' file types are allowed.",
        dictRemoveFileConfirmation: "Are you sure you want to remove this file?",
        dictMaxFilesExceeded: "The number of files uploaded exceeds the number of expected blast results.",
        retryChunks: true,
        init: function() {
          axios
            .post(
              this.options.url.slice(0,this.options.url.indexOf("drop")) + `numFiles/none`
            )
            .then((response) => {
              console.log(response.data);
              this.options.maxFiles = Number(response.data);
            })
            .catch((error) => {
              console.log(error);
            });

          axios
            .post(
              this.options.url.slice(0,this.options.url.indexOf("drop")) + `displayOutput/none`
            )
            .then((response) => {
              console.log(response.data);
              var fileNames = response.data.file_names;
              var fileSizes = response.data.file_sizes;
              for (var i = 0; i < fileNames.length; ++i) {
                this.addCustomFile(
                  // File options
                  {
                      // flag: processing is complete
                      processing: true,
                      // flag: file is accepted (for limiting maxFiles)
                      accepted: true,
                      // name of file on page
                      name: fileNames[i],
                      // image size
                      size: fileSizes[i],
                      // image type
                      type: 'application/json',
                      // flag: status upload
                      status: this.SUCCESS,
                      lastModifiedDate: "unimportant",
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

          // this.on("addedfile", function(file) {
          //   console.log(file.lastModifiedDate);
          //   file.lastModifiedDate = "unimportant";
          // });

          this.on("removedfile", function(file) {
            console.log(file);
            if (file.processing) {
            axios
              .post(
                this.options.url.slice(0,this.options.url.indexOf("drop")) + `deleteOutput/${file.name}`
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
     * @return {string} the upload URL for dropzone
     */
    getUploadUrl() {
      return `http://daniel.byu.edu:5000/phlash_api/blast/${this.$route.params.phageID}/drop/${this.numFiles}`;
    },

    autoAnnotate() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/autoAnnotate/none`
        )
        .then((response) => {
          console.log(response.data);
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Checks to see if the blast output files have been uploaded.
     * Checks to see if the blast input folder has been downloaded.
     */
    checkFiles() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/checkFiles/none`
        )
        .then((response) => {
          console.log(response.data);
          this.blastDownloaded = response.data.blast_downloaded;
          if (this.fileNames.length == this.numFiles) this.blastUploaded = true;
          if (response.data.message == "Uploaded") this.blastUploaded = true;
          console.log(this.blastUploaded);
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
          this.autoAnnotate();
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
      for (var i = 0; i < this.numFiles; ++i) {
        window.open(
          "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins",
          "_blank"
        );
      }
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
          this.fileNames = response.data.file_names;
          this.checkFiles();
        })
        .catch((error) => {
          console.log(error);
        });
    },

    uploadReminder() {
      if (!this.blastUploaded) {
        this.$bvToast.toast(`All ${this.numFiles} files have not been uploaded. If it looks like all files have been added, 
        it is possible that a duplicate file has been added. Remove the duplicate and re-upload it and all other files before continuing.`, {
          title: 'UPLOAD ALL FILES',
          autoHideDelay: 15000,
          appendToast: false
        })
      }
    },

    removeAll() {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/deleteBlastResults/none`
        )
        .then((response) => {
          if (response.data == "fail") {
            this.$bvToast.toast(`The BLAST results are currently being interpretted which may take several 
            minutes. Removal of the BLAST output files is not possible at this time.`, {
              title: 'BLAST OUTPUT FILES IN USE',
              autoHideDelay: 5000,
              appendToast: false
            })
          }
          else {
            this.blastUploaded = false;
          }
          console.log(response.data);
        })
        .catch((error) => {
          console.log(error);
        });

    }

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
