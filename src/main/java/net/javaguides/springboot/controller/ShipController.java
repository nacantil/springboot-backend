package net.javaguides.springboot.controller;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.ResponseEntity;
import org.springframework.web.bind.annotation.CrossOrigin;
import org.springframework.web.bind.annotation.DeleteMapping;
import org.springframework.web.bind.annotation.GetMapping;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.PostMapping;
import org.springframework.web.bind.annotation.PutMapping;
import org.springframework.web.bind.annotation.RequestBody;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RestController;

import net.javaguides.springboot.exception.ResourceNotFoundException;
import net.javaguides.springboot.model.Ship;
import net.javaguides.springboot.repository.ShipRepository;

import gov.nasa.worldwind.geom.LatLon;

@CrossOrigin(origins = "http://localhost:3000")
@RestController
@RequestMapping("/api/v1/")
public class ShipController {

	@Autowired
	private ShipRepository shipRepository;
	
	// get all ships
	@GetMapping("/ships")
	public List<Ship> getAllShips(){
		return shipRepository.findAll();
	}		
	
	// create ship rest api
	@PostMapping("/ships")
	public Ship createShip(@RequestBody Ship ship) {
		LatLon latLon = LatLon.fromDegrees(33.33, 11.22);
		System.out.println(latLon);
		return shipRepository.save(ship);
	}
	
	// get ship by id rest api
	@GetMapping("/ships/{id}")
	public ResponseEntity<Ship> getShipById(@PathVariable Long id) {
		Ship ship = shipRepository.findById(id)
				.orElseThrow(() -> new ResourceNotFoundException("Ship not exist with id :" + id));
		return ResponseEntity.ok(ship);
	}
	
	// update ship rest api
	
	@PutMapping("/ships/{id}")
	public ResponseEntity<Ship> updateShip(@PathVariable Long id, @RequestBody Ship shipDetails){
		Ship ship = shipRepository.findById(id)
				.orElseThrow(() -> new ResourceNotFoundException("Ship not exist with id :" + id));
		
		ship.setName(ship.getName());
		ship.setLatitude(shipDetails.getLatitude());
		ship.setLongitude(shipDetails.getLongitude());
		
		Ship updatedShip = shipRepository.save(ship);
		return ResponseEntity.ok(updatedShip);
	}
	
	// delete ship rest api
	@DeleteMapping("/ships/{id}")
	public ResponseEntity<Map<String, Boolean>> deleteShip(@PathVariable Long id){
		Ship ship = shipRepository.findById(id)
				.orElseThrow(() -> new ResourceNotFoundException("Ship not exist with id :" + id));
		
		shipRepository.delete(ship);
		Map<String, Boolean> response = new HashMap<>();
		response.put("deleted", Boolean.TRUE);
		return ResponseEntity.ok(response);
	}

	@DeleteMapping("/deleteAllShips")
	public ResponseEntity<Map<String, Boolean>> deleteAllShips(){
		List<Ship> ships = shipRepository.findAll();
		for (Ship ship : ships) {
			shipRepository.delete(ship);
		}
		Map<String, Boolean> response = new HashMap<>();
		response.put("deleted", Boolean.TRUE);
		return ResponseEntity.ok(response);
	}
}