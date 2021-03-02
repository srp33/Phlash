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
      <h1>Contact</h1>
      <div>
        <hr />
        <div
          class="alert alert-secondary"
          style="float: left; height: 9em; width: 20%"
        >
          <img
            id="step-two"
            src="/phlash/images/Piccolo.jpg"
            style="float: center; width: 100%"
          />
        </div>
        <div
          class="alert alert-secondary"
          style="float: right; height: 9em; width: 80%"
        >
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
        </div>
        <button class="btn btn-dark btn-nav" @click="goBack()">
          <strong>&#129052; Back</strong>
        </button>
        <hr />
      </div>
    </div>
  </div>
</template>

<script>
import Navbar from '../components/Navbar.vue';

export default {
  name: 'Contact',
  components: {
    Navbar,
  },

  data() {
    return {
      prevRoute: null,
    };
  },

  beforeRouteEnter(to, from, next) {
    next((vm) => {
      vm.prevRoute = from;
    });
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

    goToWebsite(site) {
      if (site === 'Home') {
        window.open('https://biology.byu.edu/piccolo-lab', '_blank');
      } else if (site === 'Contact') {
        window.open('https://biology.byu.edu/piccolo-lab/contact', '_blank');
      } else if (site === 'GitHub') {
        window.open('https://biology.byu.edu/piccolo-lab', '_blank');
      }
    },
  },
};
</script>

<style scoped>
.btn-nav {
  margin: 0.25em;
}

.alert-secondary {
  text-align: left;
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
