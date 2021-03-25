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
          <strong style="color: #4b72b9">Blue:</strong> An arrow pointing right
          indicates a gene on the direct strand.<br />
          <strong style="color: #e98b69">Orange:</strong> An arrow pointing left
          indicates a gene on the complimentary strand.<br />
        </p>
        <hr />
        <div class="nav-btns-wrapper">
          <button class="btn btn-dark btn-nav" @click="goBack()">
            <strong>&#129052; Back</strong>
          </button>
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
      <img v-bind:src="image" />
    </div>
    <div class="container">
      <div class="alert alert-secondary">
        <hr />
        <div class="nav-btns-wrapper">
          <button class="btn btn-dark btn-nav" @click="goBack()">
            <strong>&#129052; Back</strong>
          </button>
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
      pageLoading: true,
      image: null,
      prevRoute: null,
    };
  },

  beforeRouteEnter(to, from, next) {
    next((vm) => {
      vm.prevRoute = from;
    });
  },

  beforeCreate() {
    Vue.use(LoaderPlugin, {
      client_id: process.env.GOOGLE_CLIENT_ID
    });
    Vue.GoogleAuth.then(auth2 => {
      if (!auth2.isSignedIn.get()) {
        this.$router.push('/');
      }
      axios
        .get(process.env.VUE_APP_BASE_URL + `/check_user/${auth2.currentUser.get().Qs.zt}/${this.$route.params.phageID}`)
        .then((response) => {
          if (response.data === "fail") {
            this.$router.push('/');
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
      return true;
    },

    navBlast: function () {
      return true;
    },

    navAnnotations: function () {
      return true;
    },

    navGeneMap: function () {
      return true;
    },

    navSettings: function () {
      return true;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },
  },

  methods: {
    goBack() {
      console.log(this.prevRoute);
      this.$router.push(this.prevRoute);
    },

    getGraph() {
      axios
        .get(
          process.env.VUE_APP_BASE_URL +
            `/annotations/geneMap/${this.$route.params.phageID}`
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
