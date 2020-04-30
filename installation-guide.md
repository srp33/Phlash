# Phlash Installation Guide

## Back end
The back end uses [Flask](https://flask.palletsprojects.com/), a Python-based micro web framework.
### Virtual environment
First, create a virtual environment inside our project's `back-end` directory using `requirements.txt`. This file contains all of the project's required dependencies and packages.
```bash
$ cd back-end
$ python3 -m venv phlash-venv/
$ source phlash-venv/bin/activate  # activate the environment
(phlash-venv) $ pip install -r requirements.txt  # new fancy command prompt
```
> :information_source: You can replace `phlash-venv/` with a different name for your environment.

> :information_source: An *active* virtual environment is indicated by the command prompt prefixed with the name of the active environment in parentheses. 

To make sure your active environment is using the right interpreters, libraries, etc., check the value of your `PATH`. Your `PATH` should look more or less like this. 
```bash
(phlash-venv) $ echo $PATH
~/Phlash/back-end/phlash-venv/bin:/usr/local/bin:/usr/bin:/usr/sbin:/bin:/sbin
```
You can also check that your shell knows to use our project's local Python instance.
```bash
(phlash-venv) $ which python3
~/Phlash/back-end/phlash-venv/bin/python3
```

When you're done working on the project, exit the environment.
```bash
(phlash-venv) $ deactivate
$                           # old familiar command prompt
```

### Running Flask
We can now run Flask because it has been installed in our virtual environment, along with other packages. At this point, you should still be in the `back-end` directory with your environment activated. 
```bash
(phlash-venv) $ pwd
~/Phlash/back-end
```
#### Development server
Use the following commands to run the Flask server locally. Keep this running while you move on to the front end.
```bash
(phlash-venv) $ export FLASK_ENV=development  # optional
(phlash-venv) $ flask run
 * Environment: development
 * Debug mode: on
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)
```
> :information_source: If `FLASK_ENV` is set to `development`, the `flask` command will enable debug mode, and `flask run` will enable the interactive debugger and reloader.

#### Deployment options
To run in production, see [Deployment Options](https://flask.palletsprojects.com/en/1.1.x/deploying/#deployment).

## Front end
The front end uses [Vue.js](https://vuejs.org/v2/guide/), a JavaScript framework for building user interfaces.

### Running Vue.js
While running the back end in your current terminal, open another terminal and run the front end by doing the following.
```bash
$ cd front-end
$ npm install
$ npm run serve
```
> :information_source: `npm install` downloads dependencies defined in the `package.json` file and generates a `node_modules` folder with the installed modules.

You should now be able to browse to `localhost:5050` and navigate through Phlash.

