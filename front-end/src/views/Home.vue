<template>
  <div class="wrapper">
    <Navbar
      :upload="navUpload"
      :blast="navBlast"
      :annotations="navAnnotations"
      :geneMap="navGeneMap"
      :settings="navSettings"
      :phageID="navPhageID"
      :logout="loggedIn"
    />
    <div class="container">
      <div v-if="!loggedIn">
        <p style="margin: 1em">
          <img id="logo" src="/phlash/images/logohome.png" width="250" />
        </p>
        <h1>Phlash</h1>
        <div class="alert alert-secondary">
          <p style="text-align: center">
            <strong
              >A user-friendly bacteriophage genome annotation
              application.</strong
            >
          </p>
          <hr />
          <div style="text-align:center; max-height:2em;">
          <GoogleLogin class="btn btn-dark" style="margin-top: 0px; font-size: 1em;" :params="params" :onSuccess="onSuccess" :onFailure="onFailure">To get started, login with Google.</GoogleLogin>
          </div>
        </div>
      </div>
      <div v-if="loggedIn">
        <h1>Home *TEST*</h1>
        <p style="margin: 1em">
          <img id="userImage" src="/phlash/images/logohome.png" width="100" />
        </p>
        <h2>{{userName}}</h2>
        <div class="alert alert-secondary">
          <hr />
        <p v-if="phageNames.length !== 0" style="text-align: left">
          Enter an ID that contains only letters, numbers, and underscores
          below or continue with an existing ID. Each ID uniquely identifies a bacteriophage genome
          annotation.<br />
        </p>
        <p v-if="phageNames.length === 0" style="text-align: left">
          Enter an ID that contains only letters, numbers, and underscores
          below. Each ID uniquely identifies a bacteriophage genome
          annotation.<br />
        </p>
        <div
          class="input-group mb-2"
          style="float: center; width: 50%; margin: 0 auto"
        >
          <input
            class="form-control"
            style="height: 2em; font-size: 15pt"
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
              <strong>Create</strong>
            </button>
          </div>
          <hr />
        </div>
        <hr v-if="phageNames.length > 0" />
        <h4 style="text-align:center;" v-if="phageNames.length > 0">My Annotations</h4>
        <div v-if="phageNames.length !== 0"
          class="table table-responsive table-secondary"
          style="overflow-y: auto; max-height: 50em"
        >
          <table class="table table-hover" align="center">
            <thead>
              <tr>
                <th scope="col">Phage ID</th>
                <th scope="col">Creation Date</th>
                <th scope="col">Deletion Date</th>
                <th scope="col">Delete</th>
                <th scope="col">Continue Annotations</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(curr, index) in phageNames" :key="index">
                <td>{{curr}}</td>
                <td>{{phageCreationDates[index]}}</td>
                <td>{{phageDeletionDates[index]}}</td>
                <td>
                  <button
                    type="button"
                    class="btn btn-dark"
                    @click="deletePhage(index)"
                  >
                    <strong>&#128465; Delete</strong>
                  </button>
                </td>
                <td>
                  <button
                    type="button"
                    class="btn btn-dark"
                    @click="nextPage(index)"
                  >
                    <strong>Continue &#129054;</strong>
                  </button>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
        <hr v-if="phageViewNames.length > 0" />
        <h4 style="text-align:center;" v-if="phageViewNames.length > 0">Shared Annotations</h4>
        <div v-if="phageViewNames.length > 0"
          class="table table-responsive table-secondary"
          style="overflow-y: auto; max-height: 50em"
        >
          <table v-if="phageViewNames.length > 0" class="table table-hover" align="center">
            <thead>
              <tr>
                <th scope="col">Phage ID</th>
                <th scope="col">Shared By</th>
                <th scope="col">Delete</th>
                <th scope="col">View Annotations</th>
              </tr>
            </thead>
            <tbody>
              <tr v-for="(curr, index) in phageViewNames" :key="index">
                <td>{{curr}}</td>
                <td>{{phageViewEmails[index]}}</td>
                <td>
                  <button
                    type="button"
                    class="btn btn-dark"
                    @click="deleteViewPhage(index)"
                  >
                    <strong>&#128465; Delete</strong>
                  </button>
                </td>
                <td>
                  <button
                    type="button"
                    class="btn btn-dark"
                    @click="nextViewPage(index)"
                  >
                    <strong>&#128065; View</strong>
                  </button>
                </td>
              </tr>
            </tbody>
          </table>
        </div>
      </div>
        <em>Please note that all data associated with each Phage ID will be removed after 90 days.</em>
      </div>
    </div>
    <b-toast id="home-status" variant="primary" no-auto-hide>
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
import GoogleLogin from 'vue-google-login';
import { LoaderPlugin } from 'vue-google-login';
import Vue from 'vue';

export default {
  name: 'Home',
  components: {
    Navbar,
    GoogleLogin,
    LoaderPlugin,
  },

  data() {
    return {
      user: null,
      imageURL: null,
      userName: null,
      phageNames: [],
      phageDeletionDates: [],
      phageCreationDates: [],
      phageIDs: [],
      phagePages: [],
      phageViewNames: [],
      phageViewEmails: [],
      phageViewIDs: [],
      phageViewPages: [],
      loggedIn: false,
      phageID: null,
      idStatus: '',
      allFilesUploaded: false,
      dateToBeDeleted: null,
      blastComplete: false,
      params: {
        client_id: process.env.VUE_APP_API_KEY,
      },
      renderParams: {
        width: 250,
        height: 50,
        longtitle: true
      },
      statusMessage: "",
      statusTitle: "",
    };
  },

  beforeCreate() {
    Vue.use(LoaderPlugin, {
      client_id: process.env.VUE_APP_API_KEY
    });
    Vue.GoogleAuth.then(auth2 => {
      if (auth2.isSignedIn.get()) {
        this.loggedIn = true;
      }
    })
  },

  created() {
    Vue.GoogleAuth.then(auth2 => {
      if (auth2.isSignedIn.get()) {
        this.loggedIn = true;
        this.user = auth2.currentUser.get().getBasicProfile().getEmail();
        this.userName = auth2.currentUser.get().getBasicProfile().getName();
        this.imageURL = auth2.currentUser.get().getBasicProfile().getImageUrl();
        axios
          .get(process.env.VUE_APP_BASE_URL + `/get_user_data/${this.user}`)
          .then((response) => {
            this.loggedIn = true;
            if (response.data !== "empty") {
              this.phageNames = response.data.phage_id_list;
              this.phageCreationDates = response.data.phage_creation_date_list;
              this.phageDeletionDates = response.data.phage_deletion_date_list;
              this.phageIDs = response.data.id_list;
              this.phagePages = response.data.phage_pages;
              this.phageViewNames = response.data.phage_view_id_list;
              this.phageViewEmails = response.data.phage_view_email_list;
              this.phageViewIDs = response.data.id_view_list;
              this.phageViewPages = response.data.phage_view_pages;
            }
            document.getElementById('userImage').src = this.imageURL;
          })
          .catch((error) => {
            console.error(error);
          });
      }
    })
  },

  watch: {
    phageID() {
      this.phageID = this.phageID.replace(/[^a-zA-Z0-9_]/g, '');
      if (this.phageID !== null) {
        document.getElementById('start').classList.remove('disabled');
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
     * Refreshes the page when the user logs in successfully.
     * @param {object} googleUser the account that was logged in.
     */
    onSuccess(googleUser) {
      window.location.reload();
    },

    /** 
     * Handles a failed user login.
     */
    onFailure() {
      this.statusMessage = `Google sign in failed. Please try again.
                            Note that if cookies are not enabled Google sign in will not work.
                            If you are in incognito mode on your web browser, make sure to allow cookies.`;
      this.statusTitle = "LOGIN FAILED";
      this.$bvToast.show('home-status');
    },

    /**
     * Removes a phage.
     * @param {int} index the index of the phage to be removed.
     */
    deletePhage(index) {
      var cont = confirm('Are you sure you want to permanently delete this phage ID? This will remove all progress that you have made on this phage.');
      if (cont === true) {
        axios
          .delete(process.env.VUE_APP_BASE_URL + `/home/${this.user}/${this.phageIDs[index]}`)
          .then((response) => {
            this.phageNames.splice(index,1);
            this.phageCreationDates.splice(index,1);
            this.phageIDs.splice(index,1);
            this.phageDeletionDates.splice(index,1);
            this.phagePages.splice(index,1);
          })
          .catch((error) => {
            console.error(error);
          });
      }
    },

    /**
     * Removes view access for a phage.
     * @param {int} index the index of the phage to be removed.
     */
    deleteViewPhage(index) {
      var cont = confirm('Are you sure you want to remove your view access for this phage?');
      if (cont === true) {
        axios
          .delete(process.env.VUE_APP_BASE_URL + `/home/${this.user}/${this.phageViewIDs[index]}`)
          .then((response) => {
            this.phageViewNames.splice(index,1);
            this.phageViewEmails.splice(index,1);
            this.phageViewIDs.splice(index,1);
            this.phageViewPages.splice(index,1);
          })
          .catch((error) => {
            console.error(error);
          });
      }
    },

    /**
     * Navigates to the next page.
     * @param {int} index the index of the phage to be annotated.
     */
    nextPage(index) {
      this.$router.push(`/${this.phagePages[index]}/${this.phageIDs[index]}`);
    },

    nextViewPage(index) {
      this.$router.push(`/${this.phageViewPages[index]}/${this.phageViewIDs[index]}`);
    },

    /**
     * Checks for non expired phage ID.
     * Adds phage ID if non-existant.
     * @param {string} phageID the ID of the phage to be logged in or registered.
     */
    checkPhageID() {
      if (this.phageID !== '' && this.phageID !== null) {
        axios
          .post(process.env.VUE_APP_BASE_URL + `/home/${this.user}/${this.phageID}`)
          .then((response) => {
            this.idStatus = response.data.id_status;
            if (!this.idStatus.includes('ID already exists')) {
              this.phageIDs.push(response.data.id);
              this.phageNames.push(response.data.phage_id);
              this.phageDeletionDates.push(response.data.deletion_date);
              this.phageCreationDates.push(response.data.creation_date);
              this.phagePages.push(response.data.phage_page);
              this.statusMessage = `${this.idStatus}`;
              this.statusTitle = "ID CREATED";
              this.$bvToast.show('home-status');
            }
            else {
              this.statusMessage = `${this.idStatus}`;
              this.statusTitle = "ID ALREADY EXISTS";
              this.$bvToast.show('home-status');
            }
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
  margin-top: 0.7em;
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
  border-color: white;
  font-size: 1.4em;
  text-align: left;
}

.table-responsive thead th {
  position: sticky;
  top: 0;
  background: white;
  color: black;
  font-size: 1em;
}

h1 {
  margin-top: 0.7em;
}

.text-size {
  font-size: 1.2em;
}

</style>