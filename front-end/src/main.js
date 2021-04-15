import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap-vue/dist/bootstrap-vue.css';
import BootstrapVue from 'bootstrap-vue';
import Vue from 'vue';
import App from './App.vue';
import router from './router';
import GoogleSignInButton from 'vue-google-signin-button-directive';
import GSignInButton from 'vue-google-signin-button';
import GoogleLogin from 'vue-google-login';
import { LoaderPlugin } from 'vue-google-login';

Vue.use(BootstrapVue);
Vue.use(GSignInButton);
// console.log('http://daniel.byu.edu:5000/phlash_api/google_client_id.txt'.toURL().text);
console.log(process.env.VUE_APP_BASE_URL);
console.log(process.env.GOOGLE);
Vue.use(LoaderPlugin, {
  client_id: "780981769382-odbkfqn6mr1f2d9kkeaokbks7eqfrvu7.apps.googleusercontent.com"
});

Vue.GoogleAuth.then(auth2 => {
  if (!auth2.isSignedIn.get()) {
    router.push('/');
  }
  console.log(auth2.isSignedIn.get());
  console.log(auth2.currentUser.get())
})
Vue.config.productionTip = false;

new Vue({
  router,
  GoogleSignInButton,
  GSignInButton,
  GoogleLogin,
  LoaderPlugin,
  render: (h) => h(App),
}).$mount('#app');
