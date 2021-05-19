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
        This application was created by the Piccolo Lab at Brigham Young
        University. <br />
        To view the source code, report a bug, or suggest a change, click
        <a href="#" @click="goToWebsite('GitHub')" class="alert-link"
          ><i>here</i></a
        >. Pull requests are welcomed.<br />
        To contact Stephen Piccolo, click
        <a href="#" @click="goToWebsite('Contact')" class="alert-link"
          ><i>here</i></a
        >. <br />
        To view more information about the Piccolo Lab, click
        <a href="#" @click="goToWebsite('Home')" class="alert-link"
          ><i>here</i></a
        >. <br />
        <hr />
        <div class="nav-btns-wrapper">
          <button class="btn btn-dark btn-nav" @click="goBack()">
            <strong>&#129052; Back</strong>
          </button>
        </div>
        <hr />
      </div>
      <div
        class="alert alert-secondary"
        style="height: 15em; width: 100%"
      >
        <img
          src="/phlash/images/Piccolo.jpg"
          style="float: left; width: 33.33%;"
        />
        <img
          src="/phlash/images/byu.png"
          style="float: center; width: 33.33%;"
        />
        <img
          src="/phlash/images/lsb.jpeg"
          style="float: right; width: 33.33%;"
        />
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

    goToWebsite(site) {
      if (site === 'Home') {
        window.open('https://biology.byu.edu/piccolo-lab', '_blank');
      } else if (site === 'Contact') {
        window.open('https://biology.byu.edu/piccolo-lab/contact', '_blank');
      } else if (site === 'GitHub') {
        window.open('https://github.com/srp33/Phlash', '_blank');
      }
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
