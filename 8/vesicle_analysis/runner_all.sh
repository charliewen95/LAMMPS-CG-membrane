for dir in r350*/; do
  if [ -d "$dir" ]; then
    cd "$dir" || exit
    # Your commands here
    pwd
    ./runner.sh
    cd ..
  fi
done


