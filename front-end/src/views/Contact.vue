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
      <h1>About</h1>
      <div
        class="alert alert-secondary"
      >
        <hr />
        <p>This application was created by the <a href="https://biology.byu.edu/piccolo-lab" class="alert-link" target="_new">Piccolo Lab</a> at <a href="https://www.byu.edu" class="alert-link" target="_new">Brigham Young University</a>.</p>

        <p>Please use <a href="https://github.com/srp33/Phlash" class="alert-link" target="_new">GitHub</a> to view the source code, report a bug, or suggest a change. Pull requests are welcome.</p>

        <hr />
        
        <div class="nav-btns-wrapper">
          <button class="btn btn-dark btn-nav" @click="goBack()">
            <strong>&#129052; Back</strong>
          </button>
        </div>
        <hr />
      </div>
    </div>
  </div>
</template>

<script>
import Navbar from '../components/Navbar.vue';
import Vue from 'vue';
import { LoaderPlugin } from 'vue-google-login';

export default {
  name: 'Contact',
  components: {
    Navbar,
  },

  data() {
    return {
      prevRoute: null,
      loggedIn: false,
    };
  },

  beforeRouteEnter(to, from, next) {
    next((vm) => {
      vm.prevRoute = from;
    });
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

  computed: {
    navUpload: function () {
      return false;
    },

    navBlast: function () {
      return false;
    },

    navAnnotations: function () {
      return false;
    },

    navGeneMap: function () {
      return false;
    },

    navSettings: function () {
      return false;
    },

    navPhageID: function () {
      return this.$route.params.phageID;
    },
  },

  methods: {
    goBack() {
      this.$router.push(this.prevRoute);
    },
  },
};
</script>

<style scoped>
.btn-nav {
  margin: 0.25em;
}

.nav-btns-wrapper {
  text-align: center;
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
</style>
