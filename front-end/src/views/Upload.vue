<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
      :logout="true"
    />
    <div class="container">
      <h1>Upload</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p v-if="fasta">
          You have uploaded a FASTA file. Click 'Next' to continue.
        </p>
        <p v-else>
          Please upload a FASTA file
          <strong>(.fasta, .fna, .fa)</strong> containing the entire DNA
          sequence of your bacteriophage.
        </p>
        <p>
          On the next page this phage will be auto-annotated using programs such as 
          <a href="#" @click="goToWebsite('GeneMarkS')" class="alert-link"
            ><i>GeneMarkS</i></a
          >,
          <a href="#" @click="goToWebsite('Glimmer3')" class="alert-link"
            ><i>Glimmer3</i></a
          >,
          <a href="#" @click="goToWebsite('Prodigal')" class="alert-link"
            ><i>Prodigal</i></a
          >, and
          <a href="#" @click="goToWebsite('Aragorn')" class="alert-link"
            ><i>ARAGORN</i></a
          >. It is also possible to auto-annotate with <a href="#" @click="goToWebsite('Phanotate')" class="alert-link"
            ><i>PHANOTATE</i></a
          >, although this tool is not used by default.
          If you would like to edit which of these tools are used, go to the 'phlash' menu and then click on 'settings'.
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Home' }">
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Blast',
              params: { phageID: $route.params.phageID },
            }"
            :event="fasta ? 'click' : ''"
          >
            <button
              class="btn btn-dark btn-nav disabled"
              id="next-top"
              @click="uploadReminder"
            >
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>

      <div class="upload-wrapper">
        <vue-dropzone
          ref="myVueDropzone"
          id="dropzone"
          :duplicateCheck="true"
          :options="dropzoneOptions"
          :destroyDropzone="false"
        ></vue-dropzone>
      </div>
      <div class="alert alert-secondary">
        <hr />
        <div class="nav-btns-wrapper">
          <router-link :to="{ name: 'Home' }">
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Blast',
              params: { phageID: $route.params.phageID },
            }"
            :event="fasta ? 'click' : ''"
          >
            <button
              class="btn btn-dark btn-nav disabled"
              id="next-bottom"
              @click="uploadReminder"
            >
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
    <b-toast id="upload-status" variant="primary" no-auto-hide>
      <template #toast-title>
        <strong class="text-size"> {{statusTitle}} </strong>
      </template>
      <div class="text-size">{{ statusMessage }}</div>
    </b-toast>
  </div>
</template>

<script>
import axios from 'axios';
import vue2Dropzone from 'vue2-dropzone';
import 'vue2-dropzone/dist/vue2Dropzone.min.css';
import Navbar from '../components/Navbar.vue';
import { LoaderPlugin } from 'vue-google-login';
import Vue from 'vue';

export default {
  name: 'Upload',
  components: {
    VueDropzone: vue2Dropzone,
    Navbar,
  },

  data() {
    return {
      fasta: false,
      dropzoneOptions: this.setDropzone(),
      blastCompleted: false,
      interval: null,
      statusMessage: "",
      statusTitle: "",
    };
  },

  beforeCreate() {
    Vue.GoogleAuth.then(auth2 => {
      if (!auth2.isSignedIn.get()) {
        this.$router.push('/');
      }
      axios
        .get(process.env.VUE_APP_BASE_URL + `/check_user/${auth2.currentUser.get().getBasicProfile().getEmail()}/${this.$route.params.phageID}`)
        .then((response) => {
          if (response.data === "fail") {
            this.$router.push('/');
          }
          else if (response.data.view) {
            this.$router.push('/');
          }
        })
        .catch((error) => {
          console.error(error);
        });
    })
  },

  created() {
    this.checkIfFilesUploaded();
  },

  destroyed() {
    this.stopChecking();
  },

  computed: {
    navUpload: function navU() {
      return true;
    },

    navBlast: function navB() {
      return this.fasta;
    },

    navGeneMap: function navG() {
      return this.fasta;
    },

    navAnnotations: function navA() {
      return this.blastCompleted;
    },

    navSettings: function navS() {
      return true;
    },

    navPhageID: function navP() {
      return this.$route.params.phageID;
    },
  },

  watch: {
    fasta: function fastaCheck() {
      if (this.fasta) {
        document.getElementById('next-top').classList.remove('disabled');
        document.getElementById('next-bottom').classList.remove('disabled');
      } else {
        document.getElementById('next-top').classList.add('disabled');
        document.getElementById('next-bottom').classList.add('disabled');
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
        acceptedFiles: '.fasta, .fna, .fa',
        chunking: true,
        chunkSize: 1000000,
        maxFiles: 1,
        dictDefaultMessage: 'Drag FASTA file here or click to browse.',
        dictInvalidFileType:
          'Invalid file extension.',
        dictRemoveFileConfirmation:
          'Are you sure you want to remove this file? This will remove all progress that you have made on this phage.',
        dictMaxFilesExceeded: 'You can only upload one file.',
        init: function initiation() {
          const postURL = `${this.options.url.slice(
            0,
            this.options.url.indexOf('uploadFasta')
          )}display/none`;
          axios
            .post(postURL)
            .then((response) => {
              if (response.data.fasta_file !== 'Not found') {
                this.addCustomFile(
                  // File options
                  {
                    // flag: processing is complete
                    processing: true,
                    // flag: file is accepted (for limiting maxFiles)
                    accepted: true,
                    // name of file on page
                    name: response.data.fasta_file,
                    // image size
                    size: response.data.fasta_file_size,
                    // image type
                    type: '.gb',
                    // flag: status upload
                    status: this.SUCCESS,
                    lastModifiedDate: 'unimportant',
                  },
                  // Custom response for event success
                  {
                    status: 'success',
                  }
                );
              }
            })
            .catch((error) => {
              console.log(error);
            });

          this.addCustomFile = function customFile(file, response) {
            // Push file to collection
            this.files.push(file);
            // Emulate event to create interface
            this.emit('addedfile', file);
            // Add status processing to file
            this.emit('processing', file);
            // Add status success to file AND RUN EVENT success from responce
            this.emit('success', file, response, false);
            // Add status complete to file
            this.emit('complete', file);
          };

          this.on('removedfile', function removedFile(file) {
            if (file.processing) {
              axios
                .post(
                  `${this.options.url.slice(
                    0,
                    this.options.url.indexOf('uploadFasta')
                  )}delete/${file.name}`
                )
                .then((response) => {
                  console.log(response.data);
                })
                .catch((error) => {
                  console.log(error);
                });
            }
          });
        },
      };
    },

    /**
     * Returns the upload URL for dropzone.
     */
    getUploadUrl() {
      return (
        process.env.VUE_APP_BASE_URL +
        `/upload/${this.$route.params.phageID}/uploadFasta/none`
      );
    },

    /**
     * Stops the interval loop for checkIfFilesUploaded().
     */
    stopChecking() {
      if (this.fasta) {
        clearInterval(this.interval);
      }
    },

    /**
     * Checks to see if the fasta file has been uploaded.
     */
    checkIfFilesUploaded() {
      this.interval = setInterval(() => {
        if (this.$route.params.phageID !== undefined)
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
      }, 1000);
    },

    /**
     * If the next button is clicked prematurely a reminder appears.
     */
    uploadReminder() {
      if (!this.fasta) {
        this.statusMessage = `You must upload a FASTA file to continue. 
                              If you have uploaded a file and still cannot continue it is because the FASTA file is not in the correct FASTA format 
                              or your file contains a nucleic acid code other than 'A', 'C', 'T', or 'G'.`;
        this.statusTitle = "UPLOAD FASTA FILE";
        this.$bvToast.show('upload-status');
      }
      else {

      }
    },

    /**
     * Links to an external website.
     * @param {string} site the website to be redirected to.
     */
    goToWebsite(site) {
      if (site === 'GeneMarkS') {
        window.open(
          'https://academic.oup.com/nar/article/29/12/2607/1034721?login=true',
          '_blank'
        );
      } else if (site === 'Glimmer3') {
        window.open('http://ccb.jhu.edu/papers/glimmer3.pdf', '_blank');
      } else if (site === 'Aragorn') {
        window.open(
          'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC373265/',
          '_blank'
        );
      } else if (site === 'Phanotate') {
        window.open(
          'https://academic.oup.com/bioinformatics/article/35/22/4537/5480131',
          '_blank'
        );
      } else if (site === 'Prodigal') {
        window.open(
          'https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119',
          '_blank'
        );
      }
    },
  },
};
</script>

<style scoped>
.wrapper {
  margin: 0;
}

h1 {
  margin-top: 0.7em;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 0.3em;
}

.upload-wrapper {
  margin: 1em auto;
  margin-bottom: 1em;
}

.alert-secondary {
  background-color: white;
  border-color: white;
  font-size: 1.4em;
  text-align: left;
}

.btn-dark {
  font-size: 15pt;
}

.text-size {
  font-size: 1.2em;
}
</style>
