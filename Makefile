
all:
	bundle exec jekyll serve 


production:
	JEKYLL_ENV=production bundle exec jekyll build -d public



