<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :dnamaster="navDNAMaster"
      :blast="navBlast"
      :annotations="navAnnotations"
    />
    <div class="container">
      <!--<h1><strong>Phlash</strong></h1>-->
      <p>
        <img id="logo" src="/phlash/images/Logo.png" width="250" />
      </p>
      <div class="alert alert-primary" align="left">
        <p>
          Welcome to <strong><i>Phlash</i></strong
          >!
        </p>
        <p>
          Enter an ID for your bacteriphage below to get started.
          <em>Please note that your ID will be removed after 90 days.</em>
        </p>
        <div class="input-group mb-3">
          <input
            class="form-control"
            type="text"
            v-model="phageID"
            v-on:keyup.enter="checkPhageID(phageID)"
            placeholder="Enter a unique bacteriophage ID containing letters, numbers, and underscores only"
            aria-label="Enter a unique bacteriophage ID"
            aria-describedby="basic-addon2"
          />
          <div class="input-group-append">
            <button
              class="btn btn-dark btn-sm"
              type="button"
              @click="checkPhageID(id)"
            >
              <strong>Enter</strong>
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
            :to="{ name: 'DNAMaster', params: { phageID: phageID } }"
            v-if="idStatus.includes('ID already exists') && allFilesUploaded"
          >
            <button class="btn btn-light">
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
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Upload', params: { phageID: phageID } }"
            v-if="
              idStatus.includes('ID created') ||
              (idStatus.includes('ID already exists') && !allFilesUploaded)
            "
          >
            <button class="btn btn-light">
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
      phageID: "",
      idStatus: "",
      allFilesUploaded: false,
      dateToBeDeleted: null,
    };
  },
  watch: {
    phageID() {
      this.phageID = this.phageID.replace(/[^a-zA-Z0-9_]/g, "");
    },
  },
  computed: {
    navUpload: function () {
      if (this.phageID !== "") return true;
      else return false;
    },
    navDNAMaster: function () {
      if (this.phageID !== "" && this.allFilesUploaded) return true;
      return false;
    },
    navBlast: function () {
      return false;
    },
    navAnnotations: function () {
      return false;
    },
  },
  methods: {

    checkPhageID(phageID) {
      axios
        .post(process.env.VUE_APP_BASE_URL + `/home/${phageID}`)
        .then((response) => {
          this.allFilesUploaded = response.data.uploaded_all_files;
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