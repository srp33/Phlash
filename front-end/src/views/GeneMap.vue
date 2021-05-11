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
      <h1>Genome Map</h1>
      <div class="alert alert-secondary">
        <hr />
        <p><strong>Instructions</strong></p>
        <p>
          This is a visual representation of all of the genes in their current
          state.<br />
          Scroll right to see the rest of the genome.
        </p>
        <hr />
        <p><strong>Key</strong></p>
        <p>
          <strong style="color: #1b9e77">Green:</strong> An arrow pointing right
          indicates a gene on the direct strand.<br />
          <strong style="color: #d95f02">Orange:</strong> An arrow pointing left
          indicates a gene on the complimentary strand.<br />
          <strong style="color: #7570b3">Purple:</strong> An arrow pointing right 
          indicates a tRNA on the direct strand and an arrow pointing left
          indicates a gene on the complimentary strand.<br />
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'GenBank',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav" id="next-top">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
    <div style="overflow-x:auto; width=100%;">
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
      ></loading>
      <loading
        :active.sync="pageLoading"
        :is-full-page="true"
        :height="100"
        :width="100"
        >Creating Genome Map...</loading
      >
      <img v-bind:src="image" />
    </div>
    <div class="container">
      <div class="alert alert-secondary">
        <hr />
        <div class="nav-btns-wrapper">
          <router-link
            :to="{ name: 'Annotations', params: { phageID: $route.params.phageID } }"
          >
            <button class="btn btn-dark btn-nav">
              <strong>&#129052; Back</strong>
            </button>
          </router-link>
          <router-link
            :to="{
              name: 'GenBank',
              params: { phageID: $route.params.phageID },
            }"
          >
            <button class="btn btn-dark btn-nav" id="next-top">
              <strong>Next &#129054;</strong>
            </button>
          </router-link>
        </div>
        <hr />
      </div>
    </div>
  </div>
</template>

<script>
import axios from 'axios';
import Loading from 'vue-loading-overlay';
import 'vue-loading-overlay/dist/vue-loading.css';
import Navbar from '../components/Navbar.vue';
import { LoaderPlugin } from 'vue-google-login';
import Vue from 'vue';

export default {
  name: 'GeneMap',
  components: {
    Loading,
    Navbar,
  },

  data() {
    return {
      viewOnly: false,
      pageLoading: true,
      image: null,
      prevRoute: null,
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
            this.viewOnly = true
          }
        })
        .catch((error) => {
          console.error(error);
        });
    })
  },

  created() {
    this.getGraph();
  },

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

    getGraph() {
      console.log(process.env.VUE_APP_BASE_URL);
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/${this.$route.params.phageID}/geneMap`
        )
        .then((response) => {
          this.image =
            'data:image/png;base64, ' +
            response.data.image.slice(2, response.data.image.length - 1);
          this.pageLoading = false;
        })
        .catch((error) => {
          console.log(error);
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

.alert-secondary {
  background-color: white;
  border-color: white;
  font-size: 1.4em;
  text-align: left;
}

.btn-dark {
  font-size: 15pt;
}

h1 {
  margin-top: 0.7em;
}
</style>
