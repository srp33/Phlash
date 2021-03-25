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
          <GoogleLogin class="btn btn-dark" style="margin-top: 0px; font-size: 1em;" :params="params" :onSuccess="onSuccess" :onFailure="onFailure">To get started, continue with Google.</GoogleLogin>
          </div>
        </div>
      </div>
      <div v-if="loggedIn">
        <h1>Home</h1>
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
        </div>
        <!-- <hr /> -->
        <!-- <p class="id-status" v-if="idStatus !== ''">
          {{ idStatus }}
        </p>
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
        </div> -->
        <!-- <h2 v-if="phageNames.length !== 0" style="text-align:center;">Existing Annotations</h2> -->
        <!-- <h3 v-if="phageNames.length === 0">You have not started any annotations. Get started by creating a phage ID above.</h3> -->
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
                <th scope="col">Delete Phage</th>
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
                    @click="deletePhage(curr, index)"
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
      </div>
        <em>Please note that all data associated with each Phage ID will be removed after 90 days.</em>
      </div>
    </div>
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
      clientID: "780981769382-odbkfqn6mr1f2d9kkeaokbks7eqfrvu7.apps.googleusercontent.com",
      loggedIn: false,
      phageID: null,
      idStatus: '',
      allFilesUploaded: false,
      dateToBeDeleted: null,
      blastComplete: false,
      params: {
        client_id: process.env.GOOGLE_CLIENT_ID,
      },
      renderParams: {
        width: 250,
        height: 50,
        longtitle: true
      },
    };
  },

  beforeCreate() {
    Vue.use(LoaderPlugin, {
      client_id: process.env.GOOGLE_CLIENT_ID
    });
    Vue.GoogleAuth.then(auth2 => {
      if (auth2.isSignedIn.get()) {
        this.loggedIn = true;
      }
    })
  },

  created() {
    Vue.use(LoaderPlugin, {
      client_id: process.env.GOOGLE_CLIENT_ID
    });
    Vue.GoogleAuth.then(auth2 => {
      if (auth2.isSignedIn.get()) {
        this.loggedIn = true;
        console.log(this.loggedIn);
        this.user = auth2.currentUser.get().Qs.zt;
        this.userName = auth2.currentUser.get().Qs.Te;
        this.imageURL = auth2.currentUser.get().Qs.getImageUrl();
        console.log(this.imageURL);
        axios
          .get(process.env.VUE_APP_BASE_URL + `/get_user_data/${auth2.currentUser.get().Qs.zt}`)
          .then((response) => {
            console.log(response.data);
            this.loggedIn = true;
            console.log(this.loggedIn);
            if (response.data !== "empty") {
              this.phageNames = response.data.phage_id_list;
              this.phageCreationDates = response.data.phage_creation_date_list;
              this.phageDeletionDates = response.data.phage_deletion_date_list;
              this.phageIDs = response.data.id_list;
              this.phagePages = response.data.phage_pages;
            }
            document.getElementById('userImage').src = auth2.currentUser.get().Qs.getImageUrl();
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
      this.$bvToast.toast(
        `Google sign in failed. Please try again.`,
        {
          variant: 'primary',
          title: 'Login Failed',
          autoHideDelay: 15000,
          appendToast: false,
        }
      );
    },

    /**
     * Removes a phage.
     * @param {string} phageName the name of the phage to be removed.
     * @param {int} index the index of the phage to be removed.
     */
    deletePhage(phageName, index) {
      var cont = confirm('Are you sure you want to permanently delete this phage ID? This will remove all progress that you have made on this phage.');
      if (cont === true) {
        axios
          .delete(process.env.VUE_APP_BASE_URL + `/home/${this.user}/${phageName}`)
          .then((response) => {
            console.log(response.data);
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
     * Navigates to the next page.
     * @param {int} index the index of the phage to be annotated.
     */
    nextPage(index) {
      this.$router.push(`/${this.phagePages[index]}/${this.phageIDs[index]}`);
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
            console.log(response.data);
            this.idStatus = response.data.id_status;
            if (!this.idStatus.includes('ID already exists')) {
              this.phageIDs.push(response.data.id);
              this.phageNames.push(response.data.phage_id);
              this.phageDeletionDates.push(response.data.deletion_date);
              this.phageCreationDates.push(response.data.creation_date);
              this.phagePages.push(response.data.phage_page);
              this.$bvToast.toast(
                `${this.idStatus}`,
                {
                  variant: 'primary',
                  title: 'ID Created',
                  autoHideDelay: 5000,
                  appendToast: false,
                }
              );
            }
            else {
              this.$bvToast.toast(
                `${this.idStatus}`,
                {
                  variant: 'primary',
                  title: 'ID Already Exists',
                  autoHideDelay: 5000,
                  appendToast: false,
                }
              );
            }
            // const monthNames = [
            //   'January',
            //   'February',
            //   'March',
            //   'April',
            //   'May',
            //   'June',
            //   'July',
            //   'August',
            //   'September',
            //   'October',
            //   'November',
            //   'December',
            // ];
            // let date = new Date(response.data.delete_time);
            // this.dateToBeDeleted = `${
            //   monthNames[date.getUTCMonth()]
            // } ${date.getUTCDate()}, ${date.getUTCFullYear()}`;
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

/* .g-signin-button {
  display: inline-block;
  padding: 4px 8px;
  border-radius: 3px;
  background-color: #3c82f7;
  color: #fff;
  box-shadow: 0 3px 0 #0f69ff;
} */

/* .google-signin-button {
  color: white;
  background-color: red;
  height: 50px;
  font-size: 16px;
  border-radius: 10px;
  padding: 10px 20px 25px 20px;
  box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19);
} */

</style>