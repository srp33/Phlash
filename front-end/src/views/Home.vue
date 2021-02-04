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
      <p>
        <img id="logo" src="/phlash/images/Logo.png" width="250" />
      </p>
      <div class="alert alert-primary" align="left">
        <p>
          Welcome to <strong><i>Phlash</i></strong
          >!
        </p>
        <p>
          Enter a unique bacteriophage ID that only contains letters, numbers, and underscores below.<br>
          <em>Please note that all data associated with this ID will be removed after 90 days.</em>
        </p>
        <div class="input-group mb-3">
          <input
            class="form-control"
            type="text"
            v-model="phageID"
            v-on:keyup.enter="checkPhageID()"
            placeholder="Enter ID here."
            aria-label="Enter a unique bacteriophage ID"
            aria-describedby="basic-addon2"
          />
          <div class="input-group-append">
            <button
              class="btn btn-dark btn-sm"
              type="button"
              @click="checkPhageID()"
            >
              <strong>Start</strong>
            </button>
          </div>
        </div>
        <p class="id-status" v-if="idStatus !== ''">
          {{ idStatus }}
        </p>
        <div class="alert alert-warning" v-if="idStatus !== ''">
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
          You have until <strong>{{ dateToBeDeleted }}</strong> to complete
          annotations for this phage.
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Blast', params: { phageID: phageID } }"
            v-if="idStatus.includes('ID already exists') && allFilesUploaded && !blastComplete"
          >
            <button class="btn btn-light">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: phageID } }"
            v-if="idStatus.includes('ID already exists') && allFilesUploaded && blastComplete"
          >
            <button class="btn btn-light">
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
            <button class="btn btn-light">
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
    },

  },
  computed: {

    navUpload: function () {
      if (this.phageID !== null) return true;
      else return false;
    },

    // navDNAMaster: function () {
    //   if (this.phageID !== "" && this.allFilesUploaded) return true;
    //   return false;
    // },

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
      axios
        .post(process.env.VUE_APP_BASE_URL + `/home/${this.phageID}`)
        .then((response) => {
          console.log(response.data)
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
    },
    
  },
};
</script>

<style scoped>
.wrapper {
  margin: 0;
}

h1 {
  margin: 40px auto;
}

.id-status {
  margin-top: 10px;
  font-style: italic;
}

.nav-btns-wrapper {
  text-align: center;
}

.bi-arrow-right {
  margin-right: 0px;
  margin-left: 5px;
}
</style>