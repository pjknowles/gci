###################
#
### Pulls each dependency and updates hashes
#
###################

for dep in random_gen 
do
    echo "Pulling ${dep}"
    cd $dep
    # do not pull if in detached HEAD state
    if branch=$(git symbolic-ref --short -q HEAD)
    then
        git pull
    else
        echo "Deatached HEAD state. No pull."
    fi
    cd ../
done
. update_SHA1.sh
