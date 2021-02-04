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
      <h1>Upload</h1>
      <div class="alert alert-primary">
        <p><strong>Instructions</strong></p>
        <p v-if="fasta">You have uploaded a FASTA file. Click 'Next' to continue.</p>
        <p v-else>Please upload a FASTA file <strong>(.fasta, .fna)</strong> containing the entire DNA sequence of your bacteriophage.</p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Home' }">
            <button class="btn btn-light btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Blast',
              params: { phageID: $route.params.phageID },
            }"
            :event="fasta? 'click' : ''"
          >
            <button class="btn btn-light btn-nav disabled" id="next-top" @mouseenter="checkIfFilesUploaded">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
      </div>

      <div class="upload-wrapper fasta">
        <vue-dropzone ref="myVueDropzone" id="dropzone" :duplicateCheck="true" :options="dropzoneOptions" :destroyDropzone="false"></vue-dropzone>
      </div>
      <div class="alert alert-primary">
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Home' }">
            <button class="btn btn-light btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Blast',
              params: { phageID: $route.params.phageID },
            }"
            :event="fasta? 'click' : ''"
          >
            <button class="btn btn-light btn-nav disabled" id="next-bottom" @mouseenter="checkIfFilesUploaded">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
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

    navBlast: function () {
      return this.fasta
    },

    navGeneMap: function () {
      return this.fasta
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
      if (this.fasta) {
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
        acceptedFiles: ".fasta, .fna",
        chunking: false,
        maxFiles: 1,
        dictDefaultMessage: "Drag FASTA file here or click to browse.",
        dictInvalidFileType: "Only '.fasta' or '.fna' file types are allowed.",
        dictRemoveFileConfirmation: "Are you sure you want to remove this file? This will remove all progress that you have made on this phage.",
        dictMaxFilesExceeded: "You can only upload one file.",
        init: function() {

          axios
            .post(
              this.options.url.slice(0,this.options.url.indexOf("uploadFasta")) + `display/none`
            )
            .then((response) => {
              console.log(response.data);
              var fileName = response.data.fasta_file;
              var fileSize = response.data.fasta_file_size;
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
                this.options.url.slice(0,this.options.url.indexOf("uploadFasta")) + `delete/${file.name}`
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
      return `http://daniel.byu.edu:5000/phlash_api/upload/${this.$route.params.phageID}/uploadFasta/none`;
    },

    /**
     * Checks to see if the fasta file has been uploaded.
     */
    checkIfFilesUploaded() {
      console.log(this.$route.params.phageID);
      if (this.$route.params.phageID != undefined)
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/check_upload/${this.$route.params.phageID}`
        )
        .then((response) => {
          this.blastCompleted = response.data.blast_completed;
          this.fasta = response.data.fasta;
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

.alert-warning {
  margin-top: 20px;
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

.upload-btn-wrapper input[type="file"] {
  font-size: 100px;
  position: absolute;
  left: 0;
  top: 0;
  opacity: 0;
}

.fasta {
  margin-bottom: 15px;
}

</style>
