<template>
  <div class="wrapper">
    <Navbar :upload="navUpload" :dnamaster="navDNAMaster" :blast="navBlast" :annotations="navAnnotations" />
    <div class="container">
      <h1><strong>Phlash</strong></h1>
      <div class="alert alert-primary" align="left">
        <p>
          Welcome to <strong><i>Phlash</i>
          </strong>! Enter an ID for your bacteriphage below to get started.
          <em>Please note that your ID will be removed after 90 days.</em>
        </p>
        <div class="input-group mb-3">
          <input
            class="form-control"
            type="text"
            v-model="phageID"
            v-on:keyup.enter="checkPhageID(phageID)"
            placeholder="Enter a unique ID"
            aria-label="Enter a unique ID"
            aria-describedby="basic-addon2"
          />
          <div class="input-group-append">
            <button class="btn btn-dark btn-sm" type="button" @click="checkPhageID(id)">
              <strong>Enter</strong>
            </button>
          </div>
        </div>
        <p class="id-status">{{ idStatus }}</p>
        <router-link :to="{ name: 'DNAMaster', params: {phageID: phageID} }"
          v-show="allFilesUploaded">
          <button class="btn btn-light">
            <strong>Next</strong>
          </button>
        </router-link>
        <router-link :to="{ name: 'Upload', params: {phageID: phageID} }"
          v-show="!allFilesUploaded">
          <button class="btn btn-light">
            <strong>Next</strong>
          </button>
        </router-link>
      </div>
    </div>
  </div>
</template>

<script>
import axios from "axios";
import Navbar from "../components/Navbar.vue"

export default {
  name: "Home",
  components: {
    Navbar,
  },
  data() {
    return {
      phageID: "",
      idStatus: "",
      allFilesUploaded: false
    };
  },
  computed: {
    navUpload: function() {
      if (this.phageID !== "") return true;
      else return false;
    },
    navDNAMaster: function() {
      return false;
    },
    navBlast: function() {
      return false;
    },
    navAnnotations: function() {
      return false;
    },
  },
  methods: {
    checkPhageID(phageID) {
      axios.post(`http://localhost:5000/api/home/${phageID}`)
        .then(response => {
          this.idStatus = response.data.id_status;
          this.allFilesUploaded = response.data.uploaded_all_files;
          console.log(this.allFilesUploaded)
        })
        .catch(error => {
          console.error(error);
        });
    }
  }
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
</style>
