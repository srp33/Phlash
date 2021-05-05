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
      <h1>Blast</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p v-if="annotating">
          {{waitMessage}} You may
          not go to the previous or next page until your genome has been auto-annotated. Even so, you
          can continue to work on your phage genome.<br />
          Phlash relies on
          <a href="#" @click="goToWebsite('GeneMarkS')" class="alert-link"
            ><i>GeneMarkS</i></a
          >,
          <a href="#" @click="goToWebsite('Glimmer3')" class="alert-link"
            ><i>Glimmer3</i></a
          >,
          <a href="#" @click="goToWebsite('Aragorn')" class="alert-link"
            ><i>ARAGORN</i></a
          >, and
          <a href="#" @click="goToWebsite('Phanotate')" class="alert-link"
            ><i>PHANOTATE</i></a
          >
          to predict genes.
        </p>
        <p>
          BLAST is a powerful tool that will compare the genes in your phage's
          genome with a database of genes in other organisms. The BLAST results
          will help you determine the likelihood of a given gene call being
          accurate.
        </p>
        <p>Complete all of the steps below and then click 'Next'.</p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{
              name: 'Upload',
              params: { phageID: $route.params.phageID },
            }"
            :event="navUpload ? 'click' : ''"
          >
            <button class="btn btn-dark btn-nav disabled" id="back-top">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Annotations',
              params: { phageID: $route.params.phageID },
            }"
            :event="blastDownloaded && blastUploaded && autoAnnotated ? 'click' : ''"
          >
            <button
              class="btn btn-dark btn-nav disabled"
              id="next-top"
              @click="uploadReminder()"
            >
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
      <div class="steps" style="font-size: 1.25em">
        <loading
          :active.sync="orfLoading"
          :is-full-page="true"
          :height="100"
          :width="100"
          :loader="dots"
        ></loading>
        <loading
          :active.sync="orfLoading"
          :is-full-page="true"
          :height="100"
          :width="100"
          >Finding ORFs...</loading
        >
        <ol>
          <li class="step">
            <strong>Download and extract BLAST input files. </strong>
            Phlash will locate every open reading frame (ORF) on every frame of
            the phage's DNA sequence. Phlash will put the ORFS in
            <a href="#" @click="goToWebsite('Fasta')" class="alert-link"
              >multi-FASTA format</a
            >
            so they can be uploaded into BLAST. The more searches that BLAST has
            to make, the longer it will take. For this reason, Phlash separates
            the ORFs into multiple files. All of the FASTA files are then put
            into a zip file.
            <ol type="a">
              <li>
                <strong>Download:</strong>
                Download the zipped FASTA file(s) by clicking
                <a href="#" @click="downloadInputFiles" class="alert-link"
                  >here</a
                >.
                <p class="zipfile-tip">
                  &#128712; If you can't open your zip file, try using
                  <a href="#" @click="goToWebsite('Zip')" class="alert-link"
                    >7-Zip</a
                  >
                  (Windows only).
                </p>
                <p v-if="downloadLoading">Downloading BLAST input files...</p>
              </li>
              <li>
                <strong>Extract:</strong>
                You must extract all of the FASTA files from the zip folder. If
                you are unsure how to do this, watch one of the following
                tutorials:
                <a href="#" @click="goToWebsite('Windows')" class="alert-link"
                  >Windows OS</a
                >,
                <a href="#" @click="goToWebsite('Mac')" class="alert-link"
                  >Mac OS</a
                >,
                <a href="#" @click="goToWebsite('Chrome')" class="alert-link"
                  >Chrome OS</a
                >.
              </li>
            </ol>
          </li>
          <li class="step">
            <strong>BLAST the entire Phage DNA sequence.</strong>
            <div
              class="alert alert-primary alert-dismissible"
              v-if="clickedNCBI"
            >
              <a href="#" class="close" data-dismiss="alert" aria-label="close"
                >&times;</a
              >
              This will take several minutes. If it seems to be taking a long
              time, it is still probably working correctly. <br />
              If an error occurs refresh the page or try re-running that file.
              Oftentimes it will then work correctly.<br />
              For help troubleshooting,
              <a href="#" @click="goToWebsite('Help')" class="alert-link"
                >visit the FAQ</a
              >.
            </div>
            <ol type="a">
              <li>
                <strong>Go to the NCBI BLAST website:</strong> Click
                <a href="#" @click="goToWebsite('Blast')" class="alert-link">
                  here
                </a>
                to open up a tab for every BLAST input file. Please note that if
                pop-ups are blocked for Phlash, only one tab will open.
              </li>
              <li>
                <strong>Standard Protein BLAST:</strong>
                On the BLAST page ensure that the main heading says
                <i>'Standard Protein BLAST'</i>. If it does not, select the
                <i>'blastp'</i> tab.
              </li>
              <li>
                <strong>Upload File:</strong>
                In the <i>'Enter Query Sequence'</i> box click the
                <i>'Choose File'</i> button and upload one of the BLAST input
                files that was in the downloaded zip-folder. Ensure that you do
                not upload the same file twice, that you upload every file, and
                that you do not upload the zip folder as this will result in
                errors.
              </li>
              <li>
                <strong>Database:</strong>
                In the <i>'Choose Search Set'</i> box, under <i>'Database'</i>,
                choose <i>'Non-redundant protein sequences (nr)'</i> from the
                dropdown options box.
              </li>
              <li>
                <strong>Algorithm:</strong>
                In the <i>'Program Selection'</i> box, select
                <i>'blastp (protein-protein BLAST)'</i>
              </li>
              <li>
                <strong>BLAST:</strong>
                Click the <i>'BLAST'</i> button and wait for the search to
                finish.
              </li>
            </ol>
            <div class="alert alert-light alert-dismissible">
              <img
                id="step-two"
                src="/phlash/images/blast_step2.png"
                width="100%"
              />
            </div>
          </li>
          <li class="step">
            <strong>Download the Single-file JSON.</strong>
            In the top left table on the results page in BLAST, click on
            <i>'Download All.'</i> This will show you file formatting options
            for downloading your results. Choose
            <i>'Single-file JSON.'</i> After completing this process you will
            have a JSON that corresponds to each FASTA file that you
            uploaded.<br />
            <div class="alert alert-light alert-dismissible">
              <img
                id="step-three"
                src="/phlash/images/blast_step3.png"
                width="100%"
              />
            </div>
          </li>
          <li class="step" v-if="blastDownloaded">
            <div v-if="!blastUploaded">
              <strong
                >Upload your {{ numFiles }} single-file JSON BLAST
                results.</strong
              >
              You may upload all files at once or one at a time. Upload speed
              varies depending on your internet connection. If a file appears to
              get stuck while uploading simply cancel the upload and re-upload
              that file. Note that if you attempt to upload a duplicate file it
              will not upload. If you rename the files, do not include any
              special characters, including whitespace, as this may result in
              errors. <div style="color:red">Please, DO NOT leave this page while the files are uploading 
              as this may result in errors.</div>
              <vue-dropzone
                ref="myVueDropzone"
                id="dropzone"
                :duplicateCheck="true"
                :options="dropzoneOptions"
                :destroyDropzone="false"
              ></vue-dropzone>
            </div>
            <div v-if="blastUploaded">
              <strong>All BLAST files have been uploaded.</strong>
              If you wish to remove and reupload files, click
              <a href="#" @click="removeAll()" class="alert-link">here</a>.
            </div>
          </li>
        </ol>
      </div>
      <div class="alert alert-secondary">
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{
              name: 'Upload',
              params: { phageID: $route.params.phageID },
            }"
            :event="navUpload ? 'click' : ''"
          >
            <button class="btn btn-dark btn-nav disabled" id="back-bottom">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'Annotations',
              params: { phageID: $route.params.phageID },
            }"
            :event="blastDownloaded && blastUploaded && autoAnnotated ? 'click' : ''"
          >
            <button
              class="btn btn-dark btn-nav disabled"
              id="next-bottom"
              @click="uploadReminder()"
            >
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
    <b-toast id="blast-status" variant="primary" no-auto-hide>
      <template #toast-title>
        <strong class="text-size"> {{statusTitle}} </strong>
      </template>
      <div class="text-size">{{ statusMessage }}</div>
    </b-toast>
  </div>
</template>

<script>
import axios from 'axios';
import Navbar from '../components/Navbar.vue';
import Loading from 'vue-loading-overlay';
import 'vue-loading-overlay/dist/vue-loading.css';
import vue2Dropzone from 'vue2-dropzone';
import 'vue2-dropzone/dist/vue2Dropzone.min.css';
import { LoaderPlugin } from 'vue-google-login';
import Vue from 'vue';

export default {
  name: 'Blast',
  components: {
    VueDropzone: vue2Dropzone,
    Loading,
    Navbar,
  },

  data() {
    return {
      orfLoading: true,
      downloadLoading: false,
      clickedNCBI: false,
      blastDownloaded: false,
      blastUploaded: false,
      autoAnnotated: false,
      annotating: false,
      numFiles: 50,
      fileNames: [],
      badFiles: [],
      dropzoneOptions: this.setDropzone(),
      interval: null,
      waitMessage: null,
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
        .get(process.env.VUE_APP_BASE_URL + `/check_user/${auth2.currentUser.get().ft.Qt}/${this.$route.params.phageID}`)
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
    this.checkFiles();
  },

  destroyed() {
    this.stopChecking();
  },

  

  computed: {
    navUpload: function () {
      if (this.autoAnnotated && !this.annotating) return true;
      else return false;
    },

    navGeneMap: function () {
      return this.autoAnnotated;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      if (this.blastDownloaded && this.blastUploaded && this.autoAnnotated) return true;
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
      if (this.blastDownloaded && this.blastUploaded && this.autoAnnotated) {
        document.getElementById('next-top').classList.remove('disabled');
        document.getElementById('next-bottom').classList.remove('disabled');
      } else {
        document.getElementById('next-top').classList.add('disabled');
        document.getElementById('next-bottom').classList.add('disabled');
      }
    },

    blastUploaded: function () {
      if (this.blastDownloaded && this.blastUploaded && this.autoAnnotated) {
        document.getElementById('next-top').classList.remove('disabled');
        document.getElementById('next-bottom').classList.remove('disabled');
      } else {
        document.getElementById('next-top').classList.add('disabled');
        document.getElementById('next-bottom').classList.add('disabled');
      }
    },

    autoAnnotated: function () {
      if (this.autoAnnotated && !this.annotating) {
        document.getElementById('back-top').classList.remove('disabled');
        document.getElementById('back-bottom').classList.remove('disabled');
      } else {
        document.getElementById('back-top').classList.add('disabled');
        document.getElementById('back-bottom').classList.add('disabled');
      }
      if (this.blastDownloaded && this.blastUploaded && this.autoAnnotated) {
        document.getElementById('next-top').classList.remove('disabled');
        document.getElementById('next-bottom').classList.remove('disabled');
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
        chunkSize: 1000000,
        dictDefaultMessage: 'Drag files here or click to browse.',
        dictInvalidFileType: 'Only ".json" file types are allowed.',
        dictRemoveFileConfirmation:
          'Are you sure you want to remove this file?',
        dictMaxFilesExceeded:
          'The number of files uploaded exceeds the number of expected blast results.',
        retryChunks: true,
        init: function () {
          axios
            .post(
              this.options.url.slice(0, this.options.url.indexOf('drop')) +
                `numFiles/none`
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
              this.options.url.slice(0, this.options.url.indexOf('drop')) +
                `displayOutput/refresh`
            )
            .then((response) => {
              console.log(response.data);
              var fileNames = response.data.file_names;
              // fileNames = fileNames.concat(response.data.bad_files);
              var fileSizes = response.data.file_sizes;
              var fileMods = response.data.file_mods;
              for (var i = 0; i < fileNames.length; i += 1) {
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
                    size: Number(fileSizes[i]),
                    // image type
                    type: 'application/json',
                    // flag: status upload
                    status: this.SUCCESS,
                    // last modification date
                    lastModifiedDate: fileMods[i],
                  },
                  // Custom response for event success
                  {
                    status: 'success',
                  }
                );
              }
              console.log(this.files);
            })
            .catch((error) => {
              console.log(error);
            });

          this.addCustomFile = function (file, response) {
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

          this.on('addedfile', function (file) {
            console.log(file.lastModifiedDate);
            console.log(file.size);
            axios
              .post(
                this.options.url.slice(0, this.options.url.indexOf('drop')) +
                  `${file.lastModifiedDate}/${file.name}${file.size}`
              )
              .then((response) => {
                console.log(response.data);
              })
              .catch((error) => {
                console.log(error);
              });
          });

          this.on('removedfile', function (file) {
            console.log(file);
            if (file.processing) {
              axios
                .post(
                  this.options.url.slice(0, this.options.url.indexOf('drop')) +
                    `deleteOutput/${file.name}`
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
     * @return {string} the upload URL for dropzone.
     */
    getUploadUrl() {
      return (
        process.env.VUE_APP_BASE_URL +
        `/blast/${this.$route.params.phageID}/drop/${this.numFiles}`
      );
    },

    /**
     * Calls the auto-annotation method on back-end.
     */
    autoAnnotate() {
      this.annotating = true;
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
     * Stops the checkFiles interval loop.
     */
    stopChecking() {
      clearInterval(this.interval);
    },

    /**
     * Checks to see if the blast output files have been uploaded.
     * Checks to see if the blast input folder has been downloaded.
     */
    checkFiles() {
      console.log('checkFiles');
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/checkFiles/none`
        )
        .then((response) => {
          console.log(response.data);
          this.blastDownloaded = response.data.blast_downloaded;
          if (!this.blastDownloaded && !response.data.blast_input_in_progress) {
            this.createInputFiles();
          } else if (!response.data.blast_input_in_progress) {
            this.orfLoading = false;
          }
          if (response.data.uploaded) {
            this.blastUploaded = true;
          }
          if (response.data.annotated) {
            this.autoAnnotated = true;
          }
          if (!this.autoAnnotated && !response.data.annotation_in_progress) {
            this.autoAnnotate();
          }
          console.log(this.blastUploaded);
          this.setNumFiles();
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Calls the method on back-end that creates the input files needed for BLAST.
     * Updates how many files the user needs to upload.
     */
    createInputFiles() {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/createInput/none`
        )
        .then((response) => {
          console.log(response.data);
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Downloads the BLAST input files.
     */
    downloadInputFiles() {
      this.downloadLoading = true;
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/downloadInput/none`,
          FormData,
          {
            responseType: 'blob',
          }
        )
        .then((response) => {
          console.log(response.data);
          this.downloadLoading = false;
          let file_data = response.data;
          const blob = new Blob([file_data], { type: 'application/zip' });
          let link = document.createElement('a');
          link.href = window.URL.createObjectURL(blob);
          link.download = `${this.$route.params.phageID}_blast.zip`;
          link.click();
        })
        .catch((error) => {
          console.log(error);
        });
    },

    /**
     * Sets the number of blast output files that the user must upload.
     */
    setNumFiles() {
      console.log('setNumFiles');
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/numFiles/none`
        )
        .then((response) => {
          if (response.data !== 'None') {
            this.numFiles = Number(response.data);
          }
          this.displayOutputFiles();
        });
    },

    /**
     * Links to an external website.
     * @param {string} site the website to be redirected to.
     */
    goToWebsite(site) {
      if (site === 'Blast') {
        for (var i = 0; i < this.numFiles; i += 1) {
          window.open(
            'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins',
            '_blank'
          );
        }
        this.clickedNCBI = true;
      } else if (site === 'Help') {
        window.open(
          'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=Blastdocs',
          '_blank'
        );
      } else if (site === 'GeneMarkS') {
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
      } else if (site === 'Fasta') {
        window.open(
          'http://www.metagenomics.wiki/tools/fastq/multi-fasta-format#:~:text=FASTA%20is%20a%20text-file%20format%20for%20representing%20nucleotide,A%20multi-FASTA%20file%20contains%20multiple%20FASTA%20formated%20sequences.',
          '_blank'
        );
      } else if (site === 'Zip') {
        window.open('https://www.7-zip.org/', '_blank');
      } else if (site === 'Windows') {
        window.open('https://www.youtube.com/watch?v=HLBSS3JjAh0', '_blank');
      } else if (site === 'Mac') {
        window.open('https://www.youtube.com/watch?v=O-zHTZWuakg', '_blank');
      } else if (site === 'Chrome') {
        window.open('https://www.youtube.com/watch?v=OuaXF19UFsE', '_blank');
      }
    },

    /**
     * Gets all of the names of all of the blast output files that have been uploaded.
     * Gets all of the names of blast files that failed to upload correctly.
     * Determines if the needed number of files have been uploaded.
     * Determines if auto-annotation is complete.
     */
    displayOutputFiles() {
      this.interval = setInterval(() => {
        axios
          .post(
            process.env.VUE_APP_BASE_URL +
              `/blast/${this.$route.params.phageID}/displayOutput/standard`
          )
          .then((response) => {
            console.log(response.data);
            this.fileNames = response.data.file_names;
            this.badFiles = response.data.bad_files;
            if (
              this.fileNames.length == this.numFiles &&
              this.badFiles.length === 0
            ) {
              this.blastUploaded = true;
            }
            if (response.data.in_process) {
              this.annotating = true;
              if (response.data.position === 0) {
                this.waitMessage = "Phlash is currently auto-annotating your bacteriophage genome.";
              }
              else {
                this.waitMessage = "Your bacteriophage genome is currently number " + response.data.position + " in line to be auto-annotated.";
              }
            } else if (this.annotating === true || this.autoAnnotated === false) {
              this.annotating = false;
              this.autoAnnotated = true;
              if (response.data.result == "success") {
                this.statusMessage = "Phlash has finished auto-annotating this phage's genome.";
                this.statusTitle = "FINISHED";
                this.$bvToast.show('blast-status');
              }
              else if (response.data.result != "not complete") {
                this.statusMessage = `An unknown error occured while auto-annotating. 
                                      Please reupload your FASTA file and try again. 
                                      If this problem persits please contact us by visiting the 'about' page.`;
                this.statusTitle = "ERROR";
                this.$bvToast.show('blast-status');
              }
            }
            if (response.data.blast_input_complete) {
              this.numFiles = response.data.num_files;
              this.orfLoading = false;
              this.blastDownloaded = true;
              if (response.data.num_files == "error") {
                this.statusMessage = `An unknown error occured while auto-annotating. 
                                      Please reupload your FASTA file and try again. 
                                      If this problem persits please contact us by visiting the 'about' page.`;
                this.statusTitle = "ERROR";
                this.$bvToast.show('blast-status');
              }
            }
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
      if (!this.autoAnnotated) {
        this.statusMessage = "You must wait for auto-annotation to be completed before continuing.";
        this.statusTitle = "INCOMPLETE AUTO-ANNOTATION";
        this.$bvToast.show('blast-status');
      } else if (!this.blastUploaded) {
        this.statusMessage = `All ${this.numFiles} files have not been uploaded.`;
        this.statusTitle = "UPLOAD ALL FILES";
        this.$bvToast.show('blast-status');
      } else {
        clearInterval(this.interval);
      }
    },

    /**
     * Removes all of the data associated with the uploaded BLAST files.
     */
    removeAll() {
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/blast/${this.$route.params.phageID}/deleteBlastResults/none`
        )
        .then((response) => {
          if (response.data === 'fail') {
            this.statusMessage = `The BLAST results are currently being interpreted which may take several minutes. 
                                  Removal of the BLAST output files is not possible at this time.`;
            this.statusTitle = "BLAST OUTPUT FILES IN USE";
            this.$bvToast.show('blast-status');
          } else {
            this.blastUploaded = false;
          }
          console.log(response.data);
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

h1 {
  margin-top: 0.7em;
}

.alert-primary {
  margin-top: 0.5em;
}

.steps {
  text-align: left;
}

.step {
  margin-bottom: 0.5em;
}

.zipfile-tip {
  font-size: 0.9em;
  margin-bottom: 0;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 0.25em;
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
