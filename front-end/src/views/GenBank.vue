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
      <loading
        :active.sync="downloadLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <h1>GenBank</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p>
          The GenBank file can be downloaded at anytime by clicking 'Download'.
          You should not submit this file to GenBank until you have reviewed
          each CDS, made necessary adjustments, and are sure that your
          annotations are accurate. Please note that any genes that were deleted
          will not be included in the GenBank file which will result in
          renumbering of the locus tags. <br />
          To see a sample GenBank record, click
          <a href="#" @click="goToWebsite('Example')" class="alert-link"
            ><i>here</i></a
          >.<br />
          For instructions on how to submit data to GenBank, click
          <a href="#" @click="goToWebsite('Instructions')" class="alert-link"
            ><i>here</i></a
          >.<br />
        </p>
        <p>
          <button class="btn btn-dark" @click="showDownloadGenbank = true">
            <strong>&#10515; Download</strong>
          </button>
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'GeneMap', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
      <b-modal
        class="text-size"
        v-model="showDownloadGenbank"
        id="download-genbank-modal"
        ref="finishedModal"
        hide-footer
      >
        <template #modal-title>
          <div class="text-size">Download GenBank File</div>
        </template>
        <p>
          You may enter any additional information that you would like added to
          the GenBank file below, or leave it blank.
        </p>
        <hr />
        <b-form @submit="onDownloadGenBank" align="left">
          <b-form-group
            label="Phage Name:"
            label-size="lg"
            label-for="phage-name"
          >
            <b-form-input
              class="form-input"
              id="phage-name"
              type="text"
              v-model="genbankAnnotations.phageName"
              required
              placeholder="Lambda"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Source:" label-size="lg" label-for="source">
            <b-form-input
              class="form-input"
              id="source"
              type="text"
              v-model="genbankAnnotations.source"
              placeholder="Escherichia Phage Lambda"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Organism:" label-size="lg" label-for="organism">
            <b-form-input
              class="form-input"
              id="organism"
              type="text"
              v-model="genbankAnnotations.organism"
              placeholder="Viruses; Duplodnaviria; Heunggongvirae;..."
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Molecular Type:"
            label-size="lg"
            label-for="molType"
          >
            <b-form-input
              class="form-input"
              id="molType"
              type="text"
              v-model="genbankAnnotations.molType"
              placeholder="genomic DNA"
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Isolation Source:"
            label-size="lg"
            label-for="isolation-source"
          >
            <b-form-input
              class="form-input"
              id="isolation-source"
              type="text"
              v-model="genbankAnnotations.isolationSource"
              placeholder="raw sewage"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Lab Host:" label-size="lg" label-for="lab-host">
            <b-form-input
              class="form-input"
              id="lab-host"
              type="text"
              v-model="genbankAnnotations.labHost"
              placeholder="Escherichia"
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Identified By:"
            label-size="lg"
            label-for="identified-by"
          >
            <b-form-input
              class="form-input"
              id="identified-by"
              type="text"
              v-model="genbankAnnotations.identifiedBy"
              placeholder="John Doe"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Authors:" label-size="lg" label-for="authors">
            <b-form-input
              class="form-input"
              id="authors"
              type="text"
              v-model="genbankAnnotations.authors"
              placeholder="Doe,J.W., Wright,S."
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Title:" label-size="lg" label-for="title">
            <b-form-input
              class="form-input"
              id="title"
              type="text"
              v-model="genbankAnnotations.title"
              placeholder="Direct Submission"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Journal:" label-size="lg" label-for="journal">
            <b-form-input
              class="form-input"
              id="journal"
              type="text"
              v-model="genbankAnnotations.journal"
              placeholder="Submitted (18-SEP-2021) Microbiology and Molecular Biology, Brigham
              Young University, Provo, UT 84602, USA"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Country:" label-size="lg" label-for="country">
            <b-form-input
              class="form-input"
              id="country"
              type="text"
              v-model="genbankAnnotations.country"
              placeholder="USA"
            ></b-form-input>
          </b-form-group>
          <b-form-group label="Notes:" label-size="lg" label-for="notes">
            <b-form-input
              class="form-input"
              id="notes"
              type="text"
              v-model="genbankAnnotations.notes"
              placeholder="complete genome"
            ></b-form-input>
          </b-form-group>
          <b-form-group
            label="Include CDS Notes:"
            label-size="lg"
            label-for="include-notes"
          >
            <b-form-checkbox
              id="include-notes"
              type="checkbox"
              v-model="genbankAnnotations.includeNotes"
            ></b-form-checkbox>
          </b-form-group>
          <hr />
          <b-button type="submit" class="mt-3" block style="margin-top: 0em">
            <strong>Download Genbank File</strong>
          </b-button>
        </b-form>
      </b-modal>
    </div>
  </div>
</template>

<script>
import Navbar from '../components/Navbar.vue';
import axios from 'axios';
import Loading from 'vue-loading-overlay';
import 'vue-loading-overlay/dist/vue-loading.css';
import { LoaderPlugin } from 'vue-google-login';
import Vue from 'vue';

export default {
  name: 'Contact',
  components: {
    Navbar,
    Loading,
  },

  data() {
    return {
      viewOnly: false,
      prevRoute: null,
      genbankAnnotations: {
        phageName: this.navPhageID,
        source: '',
        organism: '',
        isolationSource: '',
        labHost: '',
        identifiedBy: '',
        authors: '',
        title: '',
        journal: '',
        country: 'USA',
        molType: 'genomic DNA',
        notes: 'complete genome',
        includeNotes: false,
      },
      showDownloadGenbank: false,
      downloadLoading: false,
    };
  },

  beforeRouteEnter(to, from, next) {
    next((vm) => {
      vm.prevRoute = from;
    });
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
            this.viewOnly = true;
            this.genbankAnnotations.phageName = response.data.phage_id;
          }
          else {
            this.genbankAnnotations.phageName = response.data;
          }
        })
        .catch((error) => {
          console.error(error);
        });
    })
  },

  // created() {
  //   this.genbankAnnotations.phageName = this.$route.params.phageID;
  // },

  computed: {
    navUpload: function () {
      return !this.viewOnly;
    },

    navBlast: function () {
      return !this.viewOnly;
    },

    navAnnotations: function () {
      return true;
    },

    navGeneMap: function () {
      return true;
    },

    navSettings: function () {
      return !this.viewOnly;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },
  },

  methods: {

    goToWebsite(site) {
      if (site === 'Example') {
        window.open(
          'https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html',
          '_blank'
        );
      } else if (site === 'Instructions') {
        window.open('https://www.ncbi.nlm.nih.gov/genbank/submit/', '_blank');
      }
    },

    /**
     * Downloads the GenBank file.
     */
    onDownloadGenBank(evt) {
      this.downloadLoading = true;
      console.log(this.genbankAnnotations.includeNotes);
      evt.preventDefault();
      this.$refs.finishedModal.hide();
      const payload = this.genbankAnnotations;
      axios
        .post(
          process.env.VUE_APP_BASE_URL +
            `/genbank/${this.$route.params.phageID}/none`,
          payload
        )
        .then((response) => {
          let data = response.data;
          const blob = new Blob([data]);
          let link = document.createElement('a');
          link.href = window.URL.createObjectURL(blob);
          link.download = `${this.$route.params.phageID}_modified.gb`;
          link.click();
          this.downloadLoading = false;
        });
    },
  },
};
</script>

<style scoped>
.nav-btns-wrapper {
  text-align: center;
}

.btn-nav {
  margin: 0.25em;
}

h1 {
  margin-top: 0.7em;
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

.form-input {
  height: 2em; 
  font-size: 15pt;
}
</style>
