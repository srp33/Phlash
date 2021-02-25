<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
    />
    <div class="container">
      <p style="margin:1em;">
        <img id="logo" src="/phlash/images/logohome.png" width="250" />
      </p>
      <h1>Phlash</h1>
      <div class="alert alert-secondary">
        <p style="text-align: center;">
          <strong>A user-friendly bacteriophage genome annotation application.</strong>
        </p>
        <hr />
        <p style="text-align: left;">
          Enter an ID that contains only letters, numbers, and underscores
          below. This ID will uniquely identify your bacteriophage genome
          annotation.<br />
          <em
            >Please note that all data associated with this ID will be removed
            after 90 days.</em
          >
        </p>
        <div
          class="input-group mb-2"
          style="float: center; width: 50%; margin: 0 auto"
        >
          <input
            class="form-control"
            style="height: 2em; font-size: 15pt;"
            type="text"
            v-model="phageID"
            v-on:keyup.enter="checkPhageID()"
            placeholder="Phage ID"
            aria-label="Enter a unique bacteriophage ID"
            aria-describedby="basic-addon2"
          />
          <div class="input-group-append">
            <button
              id="start"
              style="font-size: 15pt"
              class="btn btn-dark btn-sm disabled"
              type="button"
              @click="checkPhageID()"
            >
              <strong>Start</strong>
            </button>
          </div>
        </div>
        <hr />
        <p class="id-status" v-if="idStatus !== ''">
          {{ idStatus }}
        </p>
        <div class="alert alert-primary" v-if="idStatus !== ''">
          &#128712; You have until <strong>{{ dateToBeDeleted }}</strong> to
          complete annotations for this phage.
        </div>
        <hr v-if="idStatus !== ''" />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: phageID } }"
            v-if="
              idStatus.includes('ID already exists') &&
              allFilesUploaded &&
              !blastComplete
            "
          >
            <button class="btn btn-dark">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: phageID } }"
            v-if="
              idStatus.includes('ID already exists') &&
              allFilesUploaded &&
              blastComplete
            "
          >
            <button class="btn btn-dark">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Upload', params: { phageID: phageID } }"
            v-if="
              idStatus.includes('ID created') ||
              (idStatus.includes('ID already exists') && !allFilesUploaded)
            "
          >
            <button class="btn btn-dark">
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
import Navbar from "../components/Navbar.vue";

export default {
  name: "Home",
  components: {
    Navbar,
  },

  data() {
    return {
      phageID: null,
      idStatus: "",
      allFilesUploaded: false,
      dateToBeDeleted: null,
      blastComplete: false,
    };
  },

  watch: {
    phageID() {
      this.phageID = this.phageID.replace(/[^a-zA-Z0-9_]/g, "");
      if (this.phageID != null) {
        document.getElementById("start").classList.remove("disabled");
      }
    },
  },
  computed: {
    navUpload: function () {
      if (this.phageID !== null) return true;
      else return false;
    },

    navBlast: function () {
      if (this.phageID !== null && this.allFilesUploaded) return true;
      return false;
    },

    navAnnotations: function () {
      return this.blastComplete;
    },

    navGeneMap: function () {
      if (this.phageID !== null && this.allFilesUploaded) return true;
      return false;
    },

    navSettings: function () {
      return false;
    },

    navPhageID: function () {
      return this.phageID;
    },
  },

  methods: {
    /**
     * Checks for non expired phage ID.
     * Adds phage ID if non-existant.
     * @param {string} phageID the ID of the phage to be logged in or registered.
     */
    checkPhageID() {
      if (this.phageID != "" && this.phageID != null) {
        axios
          .post(process.env.VUE_APP_BASE_URL + `/home/${this.phageID}`)
          .then((response) => {
            console.log(response.data);
            this.allFilesUploaded = response.data.uploaded_all_files;
            this.blastComplete = response.data.blast_complete;
            this.idStatus = response.data.id_status;
            const monthNames = [
              "January",
              "February",
              "March",
              "April",
              "May",
              "June",
              "July",
              "August",
              "September",
              "October",
              "November",
              "December",
            ];
            let date = new Date(response.data.delete_time);
            this.dateToBeDeleted = `${
              monthNames[date.getUTCMonth()]
            } ${date.getUTCDate()}, ${date.getUTCFullYear()}`;
          })
          .catch((error) => {
            console.error(error);
          });
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
  margin-top: .7em;
}

.id-status {
  margin-top: 0.5em;
  font-style: italic;
}

.nav-btns-wrapper {
  text-align: center;
}

.btn-dark {
  font-size: 15pt;
}

.alert-secondary {
  background-color: white;
  border-color:white;
  font-size: 1.40em;
  text-align: left;
}
</style>